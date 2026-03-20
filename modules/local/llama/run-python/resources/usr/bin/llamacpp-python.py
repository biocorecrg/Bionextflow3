#!/usr/bin/env python3

import argparse
import json
import re

import llama_cpp


def estimate_tokens(text):
    """Rough token estimate: ~4 chars per token"""
    return len(text) // 4


def chunk_text(text, max_tokens):
    """Split text into chunks by paragraphs/sections with hard size limit"""
    # Split by double newlines (paragraphs/sections)
    sections = re.split(r'\n\n+', text)
    
    chunks = []
    current_chunk = []
    current_length = 0
    
    for section in sections:
        section_tokens = estimate_tokens(section)
        
        # If a SINGLE section exceeds max_tokens, split it by lines
        if section_tokens > max_tokens:
            lines = section.split('\n')
            temp_chunk = []
            temp_length = 0
            
            for line in lines:
                line_tokens = estimate_tokens(line)
                
                # If even a single line is too big, split it by characters
                if line_tokens > max_tokens:
                    chars_per_chunk = max_tokens * 4  # ~4 chars per token
                    for i in range(0, len(line), chars_per_chunk):
                        chunk_piece = line[i:i+chars_per_chunk]
                        if current_chunk:
                            chunks.append('\n\n'.join(current_chunk))
                            current_chunk = []
                            current_length = 0
                        chunks.append(chunk_piece)
                    continue
                
                if temp_length + line_tokens > max_tokens:
                    if temp_chunk:
                        if current_chunk:
                            chunks.append('\n\n'.join(current_chunk))
                            current_chunk = []
                            current_length = 0
                        chunks.append('\n'.join(temp_chunk))
                    temp_chunk = [line]
                    temp_length = line_tokens
                else:
                    temp_chunk.append(line)
                    temp_length += line_tokens
            
            if temp_chunk:
                if current_chunk:
                    chunks.append('\n\n'.join(current_chunk))
                    current_chunk = []
                    current_length = 0
                chunks.append('\n'.join(temp_chunk))
            continue
        
        # Normal section handling
        if current_length + section_tokens > max_tokens:
            if current_chunk:
                chunks.append('\n\n'.join(current_chunk))
            current_chunk = [section]
            current_length = section_tokens
        else:
            current_chunk.append(section)
            current_length += section_tokens
    
    # Don't forget the last chunk
    if current_chunk:
        chunks.append('\n\n'.join(current_chunk))
    
    return chunks


def chunk_by_semantic_sections(text, max_tokens):
    """Split MultiQC-like reports by tool/module sections with size validation"""
    chunks = {}
    current_section = "General"
    current_content = []
    
    for line in text.split('\n'):
        # Detect section headers
        if re.match(r'^[A-Z][A-Za-z\s]{3,}:?\s*$', line.strip()) or \
           re.match(r'^#+\s+[A-Z]', line.strip()):
            # Save previous section
            if current_content:
                content_text = '\n'.join(current_content)
                if content_text.strip():
                    chunks[current_section] = content_text
            # Start new section
            current_section = line.strip().rstrip(':').lstrip('#').strip()
            current_content = []
        else:
            current_content.append(line)
    
    # Save last section
    if current_content:
        content_text = '\n'.join(current_content)
        if content_text.strip():
            chunks[current_section] = content_text
    
    # Validate and split oversized sections
    validated_chunks = {}
    for section_name, content in chunks.items():
        section_tokens = estimate_tokens(content)
        
        if section_tokens > max_tokens:
            print(f"  Section '{section_name}' is too large ({section_tokens} tokens). Splitting...")
            # Split this section into sub-chunks
            sub_chunks = chunk_text(content, max_tokens=max_tokens)
            for i, sub_chunk in enumerate(sub_chunks):
                validated_chunks[f"{section_name} (part {i+1})"] = sub_chunk
        else:
            validated_chunks[section_name] = content
    
    return validated_chunks


def process_with_llm(llm, prompt, temperature, max_tokens, system_prompt):
    """Helper to call LLM with a single prompt"""
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": prompt}
    ]
    
    response = llm.create_chat_completion(
        messages=messages,
        temperature=temperature,
        max_tokens=max_tokens,
    )
    
    try:
        return response["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError):
        return json.dumps(response, indent=2)


def llamacpp_python(
    messages_file=None,
    model_file=None,
    temperature=0.7,
    output="output.txt",
    context_size=4096,
    chat_format="chatml",
    seed=None,
    n_gpu_layers=0,
    max_tokens=1024,
    system_prompt="""You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields. 
    You are given data from such a report. Your task is to analyse the data, and
    give 1-2 bullet points of a very short and concise overall summary for the results.
    Don't waste words: mention only the important QC issues. If there are no issues, just say so.
    Just print one or two bullet points, nothing else.
    Please do not add any extra headers to the response.
    """,
    chunk_method="auto",
    chunk_fraction=0.5,  # Use 50% of context for each chunk by default
    input_threshold=0.75,  # Trigger chunking if input exceeds 75% of context
):
    # Read input text
    with open(messages_file, "r") as f:
        text = f.read()
    
    # Initialize model once
    print(f"Loading model: {model_file}")
    llm = llama_cpp.Llama(
        model_path=model_file, 
        chat_format=chat_format, 
        n_ctx=context_size, 
        seed=seed,
        n_gpu_layers=n_gpu_layers,
        verbose=False
    )
    print("Model loaded successfully\n")
    
    # Calculate token budgets
    estimated_tokens = estimate_tokens(text)
    system_tokens = estimate_tokens(system_prompt)
    
    # Reserve space for: system prompt + response + prompt template overhead (~300 tokens)
    overhead = system_tokens + max_tokens + 300
    
    # Maximum tokens we can use for input chunks
    max_chunk_tokens = int(context_size * chunk_fraction)
    
    # Make sure we leave room for overhead
    if max_chunk_tokens + overhead > context_size:
        max_chunk_tokens = context_size - overhead
        print(f"⚠️  Adjusted chunk size to {max_chunk_tokens} tokens to fit overhead\n")
    
    print(f"Context size: {context_size} tokens")
    print(f"Chunk budget: {max_chunk_tokens} tokens ({chunk_fraction*100:.0f}% of context)")
    print(f"Input size: ~{estimated_tokens} tokens\n")
    
    # Check if chunking is needed
    if estimated_tokens > (context_size * input_threshold):
        print(f"Input exceeds {input_threshold*100:.0f}% of context. Chunking enabled.\n")
        
        # Choose chunking method
        if chunk_method == "semantic":
            chunks_dict = chunk_by_semantic_sections(text, max_tokens=max_chunk_tokens)
            chunks = list(chunks_dict.items())
            print(f"Split into {len(chunks)} semantic sections:")
            for name, _ in chunks:
                print(f"  - {name}")
        else:  # auto
            chunks_list = chunk_text(text, max_tokens=max_chunk_tokens)
            chunks = [(f"Part {i+1}", chunk) for i, chunk in enumerate(chunks_list)]
            print(f"Split into {len(chunks)} chunks")
        
        print()
        
        # Process each chunk
        print("="*60)
        print("Processing chunks...")
        print("="*60 + "\n")
        
        summaries = []
        for i, (section_name, chunk) in enumerate(chunks):
            chunk_tokens = estimate_tokens(chunk)
            print(f"[{i+1}/{len(chunks)}] {section_name}: ~{chunk_tokens} tokens")
            
            prompt = f"""Summarize this section of a sequencing QC report. Focus on:
- Key metrics and statistics
- Any quality issues or warnings  
- Pass/fail status
- Critical findings

Section: {section_name}

{chunk}

Provide a concise summary:"""
            
            try:
                summary = process_with_llm(llm, prompt, temperature, max_tokens, system_prompt)
                summaries.append(f"### {section_name}\n{summary}")
                print(f"  ✓ Generated ({estimate_tokens(summary)} tokens)\n")
            except Exception as e:
                print(f"  ✗ Error: {e}\n")
                summaries.append(f"### {section_name}\n[Error: {str(e)}]")
        
        # Final synthesis
        print("="*60)
        print("Generating final synthesis...")
        print("="*60 + "\n")
        
        combined_summaries = "\n\n".join(summaries)
        combined_tokens = estimate_tokens(combined_summaries)
        print(f"Combined summaries: ~{combined_tokens} tokens")
        
        # Check if we need second-level reduction
        if combined_tokens > max_chunk_tokens:
            print("Combined summaries too large. Applying second-level reduction...\n")
            
            # Group into batches
            batch_size = max(len(summaries) // 3, 1)
            batches = [summaries[i:i+batch_size] for i in range(0, len(summaries), batch_size)]
            
            print(f"Batching {len(summaries)} summaries into {len(batches)} groups\n")
            
            meta_summaries = []
            for i, batch in enumerate(batches):
                print(f"[{i+1}/{len(batches)}] Synthesizing batch {i+1}...")
                batch_text = "\n\n".join(batch)
                prompt = f"Synthesize these QC summaries into a brief overview:\n\n{batch_text}"
                
                try:
                    meta_summary = process_with_llm(llm, prompt, temperature, max_tokens, system_prompt)
                    meta_summaries.append(meta_summary)
                    print(f"  ✓ Generated ({estimate_tokens(meta_summary)} tokens)\n")
                except Exception as e:
                    print(f"  ✗ Error: {e}\n")
                    meta_summaries.append(f"[Error: {str(e)}]")
            
            combined_summaries = "\n\n".join(meta_summaries)
            print(f"After reduction: ~{estimate_tokens(combined_summaries)} tokens\n")
        
        final_prompt = f"""
    You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields. 
    You are given data from such a report. Your task is to analyse the data, and
    give 1-2 bullet points of a very short and concise overall summary for the results.
    Don't waste words: mention only the important QC issues. If there are no issues, just say so.
    Just print one or two bullet points, nothing else.
    Please do not add any extra headers to the response.

    {combined_summaries}

    """
        
        print("Generating final assessment...")
        try:
            final_summary = process_with_llm(llm, final_prompt, temperature, max_tokens * 2, system_prompt)
            print(f"✓ Generated ({estimate_tokens(final_summary)} tokens)\n")
        except Exception as e:
            print(f"✗ Error: {e}\n")
            final_summary = f"[Error: {str(e)}]"
        
        # Build final report — only the overall assessment to keep output concise
        reply = final_summary
    
    else:
        # Process without chunking
        print("Input fits in context. Processing without chunking.\n")
        messages_json = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": text}
        ]
        
        try:
            response = llm.create_chat_completion(
                messages=messages_json,
                temperature=temperature,
                max_tokens=max_tokens,
            )
            reply = response["choices"][0]["message"]["content"]
        except Exception as e:
            print(f"Error: {e}")
            reply = f"Error: {str(e)}"
    
    # Write output
    with open(output, "w") as f:
        f.write(reply)
    
    print("="*60)
    print(f"✓ Output written to: {output}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Process long documents with LLM using automatic chunking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - chunks use 50% of context by default
  %(prog)s -s report.txt -m model.gguf -o summary.txt
  
  # Use 70% of context for chunks (more context per chunk)
  %(prog)s -s report.txt -m model.gguf --chunk_fraction 0.7
  
  # Semantic chunking for MultiQC reports
  %(prog)s -s llms-full.txt -m model.gguf --chunk_method semantic
  
  # Large context model - use 60% for chunks
  %(prog)s -s report.txt -m model.gguf -c 32768 --chunk_fraction 0.6
        """
    )
    
    # Required
    parser.add_argument("-s", "--messages", required=True, 
                       help="Plain text input file")
    parser.add_argument("-m", "--model", required=True, 
                       help="Path to GGUF model file")
    
    # Model parameters
    parser.add_argument("-t", "--temperature", type=float, default=0.7, 
                       help="Temperature (default: 0.7)")
    parser.add_argument("-c", "--context", type=int, default=4096, 
                       help="Context size (default: 4096)")
    parser.add_argument("-g", "--gpu_layers", type=int, default=-1, 
                       help="GPU layers, -1 for all (default: -1)")
    parser.add_argument("--max_tokens", type=int, default=1024, 
                       help="Max tokens in response (default: 1024)")
    parser.add_argument("--chat_format", default="chatml", 
                       help="Chat format (default: chatml)")
    parser.add_argument("--seed", type=int, default=None, 
                       help="Random seed")
    
    # Output
    parser.add_argument("-o", "--output", default="output.txt", 
                       help="Output file (default: output.txt)")
    
    # Chunking
    parser.add_argument("--chunk_method", choices=["auto", "semantic"], default="auto", 
                       help="Chunking method (default: auto)")
    parser.add_argument("--chunk_fraction", type=float, default=0.5,
                       help="Fraction of context to use per chunk, 0.1-0.9 (default: 0.5)")
    parser.add_argument("--input_threshold", type=float, default=0.75,
                       help="Trigger chunking if input exceeds this fraction of context (default: 0.75)")
    
    # System prompt
    parser.add_argument("--system_prompt", 
                       default="You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields.", 
                       help="System prompt")

    args = parser.parse_args()
    
    # Validate chunk_fraction
    if not 0.1 <= args.chunk_fraction <= 0.9:
        parser.error("--chunk_fraction must be between 0.1 and 0.9")

    llamacpp_python(
        messages_file=args.messages,
        model_file=args.model,
        temperature=args.temperature,
        output=args.output,
        context_size=args.context,
        chat_format=args.chat_format,
        seed=args.seed,
        n_gpu_layers=args.gpu_layers,
        max_tokens=args.max_tokens,
        system_prompt=args.system_prompt,
        chunk_method=args.chunk_method,
        chunk_fraction=args.chunk_fraction,
        input_threshold=args.input_threshold,
    )


if __name__ == "__main__":
    main()
