#!/usr/bin/env python3

import argparse
import json

import llama_cpp


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
    system_prompt="You are an expert in bioinformatics, sequencing technologies, genomics data analysis, and adjacent fields.",
):
    messages_json = []
    # Add system prompt for coherent output
    if system_prompt:
        messages_json.append({"role": "system", "content": system_prompt})
    # Read plain text file and add as user message
    if messages_file:
        with open(messages_file, "r") as f:
            text = f.read()
            messages_json.append({"role": "user", "content": text})

    llm = llama_cpp.Llama(
        model_path=model_file, 
        chat_format=chat_format, 
        n_ctx=context_size, 
        seed=seed,
        n_gpu_layers=n_gpu_layers,  # ADD THIS
        verbose=True  # ADD THIS to see GPU offloading messages
    )
    response = llm.create_chat_completion(
        messages=messages_json,
        temperature=temperature,
        max_tokens=max_tokens,
    )

    # Extract just the assistant's reply text
    try:
        reply = response["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError):
        reply = json.dumps(response, indent=2)

    with open(output, "w") as f:
        f.write(reply)


def main():
    parser = argparse.ArgumentParser(description="Submit a process with model.")
    parser.add_argument("-s", "--messages", required=True, help="Plain text prompt file")
    parser.add_argument("-m", "--model", required=True, help="Model used")
    parser.add_argument("-t", "--temperature", type=float, default=0.9, help="Temperature")
    parser.add_argument("-o", "--output", default="output.txt", help="Output text")
    parser.add_argument("-c", "--context", type=int, default=4096, help="Context size")
    parser.add_argument("--chat_format", default="chatml", help="Chat format (chatml, llama-2, llama-3, gemma...)")
    parser.add_argument("--seed", type=int, default=None, help="Defined seed")
    parser.add_argument("-g", "--gpu_layers", type=int, default=-1, help="Number of GPU layers (-1 for all)")
    parser.add_argument("--max_tokens", type=int, default=1024, help="Max tokens in response")
    parser.add_argument("--system_prompt", default="You are a helpful assistant. Answer clearly and concisely.", help="System prompt")

    args = parser.parse_args()

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
    )


if __name__ == "__main__":
    main()
