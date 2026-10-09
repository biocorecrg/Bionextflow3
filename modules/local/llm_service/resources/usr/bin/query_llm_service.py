#!/usr/bin/env python3
"""Query an LLM web service with the content of a JSON file.

Generic client: sends the JSON (as a chat-completion prompt, or as the raw
request body with --raw-json) to an LLM endpoint and saves the response as
markdown/text and raw JSON. Supports OpenAI-compatible chat completion
endpoints (vLLM, Ollama, llama-server) and custom CRG LLM endpoints.
Rendering of the response (e.g. MultiQC) is done by other tools.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request


def parse_args():
    parser = argparse.ArgumentParser(
        description="Send JSON data and prompts to an LLM service."
    )
    parser.add_argument(
        "--json-file",
        "-i",
        required=True,
        help="Path to input JSON file (e.g., bcl2fastq Stats.json).",
    )
    parser.add_argument(
        "--url",
        "-u",
        default=os.environ.get(
            "LLM_API_URL",
            os.environ.get(
                "OPENAI_BASE_URL", "http://localhost:8080/v1/chat/completions"
            ),
        ),
        help="LLM web service endpoint URL (default: $LLM_API_URL or $OPENAI_BASE_URL).",
    )
    parser.add_argument(
        "--model",
        "-m",
        default=os.environ.get("LLM_MODEL", ""),
        help="LLM model identifier (default: $LLM_MODEL or '').",
    )
    parser.add_argument(
        "--raw-json",
        action="store_true",
        help="POST only the full content of the JSON file as request body "
        "(custom API endpoints), without chat-completion wrapping.",
    )
    parser.add_argument(
        "--system-prompt",
        "-s",
        help="System prompt text or path to a system prompt file.",
    )
    parser.add_argument(
        "--user-prompt",
        "-p",
        help="Custom user prompt template. If provided, JSON content will be appended.",
    )
    parser.add_argument(
        "--output-md",
        "-o",
        default="llm_report.md",
        help="Output markdown file path for LLM response text (default: llm_report.md).",
    )
    parser.add_argument(
        "--output-json",
        "-j",
        default="llm_response.json",
        help="Output JSON file path for raw API response (default: llm_response.json).",
    )
    parser.add_argument(
        "--temperature",
        "-t",
        type=float,
        default=0.2,
        help="Sampling temperature (default: 0.2).",
    )
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=2048,
        help="Maximum tokens to generate (default: 2048).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=180,
        help="HTTP timeout in seconds (default: 180).",
    )
    return parser.parse_args()


def load_text_or_file(text_or_path: str) -> str:
    if not text_or_path:
        return ""
    if os.path.isfile(text_or_path):
        with open(text_or_path, "r", encoding="utf-8") as f:
            return f.read().strip()
    return text_or_path.strip()


def extract_content(response_data: dict) -> str:
    """Extract message content from OpenAI or alternative API response structure."""
    if "result" in response_data and isinstance(response_data["result"], str):
        return response_data["result"]
    if "choices" in response_data and len(response_data["choices"]) > 0:
        choice = response_data["choices"][0]
        if "message" in choice and "content" in choice["message"]:
            return choice["message"]["content"]
        if "text" in choice:
            return choice["text"]
    if "response" in response_data:
        return response_data["response"]
    if "content" in response_data:
        return response_data["content"]
    # Fallback to pretty-printed json if unknown format
    return json.dumps(response_data, indent=2)


def main():
    args = parse_args()

    # Read input JSON file
    try:
        with open(args.json_file, "r", encoding="utf-8") as f:
            input_json_data = json.load(f)
            formatted_json_str = json.dumps(
                input_json_data, separators=(",", ":"), ensure_ascii=False
            )
    except Exception as e:
        print(
            f"Error reading JSON file '{args.json_file}': {e}", file=sys.stderr
        )
        sys.exit(1)

    # Determine system prompt (optional)
    system_prompt = (
        load_text_or_file(args.system_prompt)
        if args.system_prompt
        else ""
    )

    # Determine user message content
    if args.user_prompt:
        user_prompt_text = load_text_or_file(args.user_prompt)
        if user_prompt_text:
            user_content = f"{user_prompt_text}\n\n{formatted_json_str}"
        else:
            user_content = formatted_json_str
    else:
        user_content = formatted_json_str

    # Normalize URL if needed
    url = args.url.strip()
    if not url.startswith("http://") and not url.startswith("https://"):
        url = f"http://{url}"

    if args.raw_json:
        # Custom API: send only the full JSON content as body
        payload = input_json_data
    else:
        # Build messages list
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": user_content})

        # Build payload
        payload = {
            "messages": messages,
            "temperature": args.temperature,
            "max_tokens": args.max_tokens,
        }
        if args.model:
            payload["model"] = args.model

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    # API key only from environment (optional, never from command line)
    api_key = os.environ.get("LLM_API_KEY") or os.environ.get("OPENAI_API_KEY")
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    req_data = json.dumps(
        payload, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    req = urllib.request.Request(
        url, data=req_data, headers=headers, method="POST"
    )

    print(
        f"[DEBUG] args.url (raw): {args.url!r}\n"
        f"[DEBUG] URL used: {url}\n"
        f"[DEBUG] raw_json mode: {args.raw_json}, body size: {len(req_data)} bytes, "
        f"Authorization header: {'yes' if 'Authorization' in headers else 'no'}",
        file=sys.stderr,
    )
    # Debug: dump request body and print the equivalent curl command
    debug_body_file = "debug_request_body.json"
    with open(debug_body_file, "wb") as f:
        f.write(req_data)
    curl_headers = " ".join(
        f"-H '{k}: {'Bearer ***' if k == 'Authorization' else v}'"
        for k, v in headers.items()
    )
    print(
        f"[DEBUG] curl equivalent:\n"
        f"curl -v -X POST '{url}' {curl_headers} "
        f"--data-binary @{debug_body_file} --max-time {args.timeout}",
        file=sys.stderr,
    )
    print(f"Calling LLM web service at: {url} (model: {args.model})...")
    try:
        with urllib.request.urlopen(req, timeout=args.timeout) as resp:
            status_code = resp.getcode()
            response_body = resp.read().decode("utf-8")
            response_json = json.loads(response_body)
    except urllib.error.HTTPError as e:
        error_msg = e.read().decode("utf-8", errors="replace")
        print(f"HTTP Error {e.code}: {e.reason}\n{error_msg}", file=sys.stderr)
        sys.exit(1)
    except urllib.error.URLError as e:
        print(f"URL Error: {e.reason}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error querying LLM service: {e}", file=sys.stderr)
        sys.exit(1)

    # Save raw response JSON
    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(response_json, f, indent=2)

    # Extract markdown summary
    content = extract_content(response_json)
    with open(args.output_md, "w", encoding="utf-8") as f:
        f.write(content.strip() + "\n")

    print(f"Successfully received response (HTTP {status_code}).")
    print(f"Report saved to: {args.output_md}")
    print(f"Raw response saved to: {args.output_json}")


if __name__ == "__main__":
    main()

