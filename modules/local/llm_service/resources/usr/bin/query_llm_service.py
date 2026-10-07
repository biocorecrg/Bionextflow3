#!/usr/bin/env python3
"""Query an LLM web service with JSON input data (e.g. demultiplexing stats).

Supports OpenAI-compatible chat completion endpoints (vLLM, Ollama,
llama-server, CRG internal LLM endpoints).
"""

import argparse
import html
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
        "--api-key",
        "-k",
        default=(
            os.environ.get("LLM_API_KEY")
            or os.environ.get("OPENAI_API_KEY")
            or ""
        ),
        help="API key / Bearer token (default: $LLM_API_KEY or $OPENAI_API_KEY).",
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
        "--output-mqc",
        default="",
        help="Output HTML file path for MultiQC custom content report.",
    )
    parser.add_argument(
        "--sample-name",
        default="Demultiplexing",
        help="Sample or run name to display in MultiQC report (default: Demultiplexing).",
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
            formatted_json_str = json.dumps(input_json_data, indent=2)
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
    if args.api_key:
        headers["Authorization"] = f"Bearer {args.api_key}"

    req_data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=req_data, headers=headers, method="POST"
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

    # Generate MultiQC custom content HTML if requested
    if args.output_mqc:
        generate_multiqc_html(content, args.output_mqc, args.sample_name, input_data=input_json_data)
        print(f"MultiQC section saved to: {args.output_mqc}")

    print(f"Successfully received response (HTTP {status_code}).")
    print(f"Report saved to: {args.output_md}")
    print(f"Raw response saved to: {args.output_json}")


def _format_cell(val, col_name: str = "", is_first: bool = False) -> str:
    col_lower = col_name.lower()
    if isinstance(val, list):
        if not val:
            return '<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; color:#64748b;">-</td>'
        count = len(val)
        items_str = ", ".join(str(x) for x in val)
        escaped_items = html.escape(items_str)
        if count <= 3:
            return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#334155;">{escaped_items}</td>'
        else:
            cell_content = (
                f'<div style="margin-bottom:4px;">'
                f'<span style="padding:2px 7px; border-radius:10px; background:#0284c7; color:#fff; font-size:11px; font-weight:600;">{count} items</span>'
                f'</div>'
                f'<div style="max-height:95px; overflow-y:auto; font-size:12px; line-height:1.45; color:#475569; word-break:break-word; background:#f8fafc; padding:6px 8px; border-radius:4px; border:1px solid #e2e8f0; font-family:monospace;">{escaped_items}</div>'
            )
            return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; min-width:220px;">{cell_content}</td>'

    val_str = str(val) if val is not None else ""
    escaped = html.escape(val_str)

    if "sample" in col_lower and "affected" not in col_lower and "barcode" not in col_lower:
        nowrap = " white-space:nowrap;" if len(val_str) < 30 else ""
        return f'<td style="padding:9px 12px; font-weight:600; color:#1e293b; border-bottom:1px solid #e2e8f0;{nowrap}">{escaped}</td>'
    elif "expected" in col_lower:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#475569; white-space:nowrap;">{escaped}</td>'
    elif any(w in col_lower for w in ["unassigned", "observed", "barcode"]) and "affected" not in col_lower:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#0284c7; font-weight:600; white-space:nowrap;">{escaped}</td>'
    elif "relation" in col_lower or "similarity" in col_lower or "issue" in col_lower:
        return f'<td style="padding:9px 12px; color:#b91c1c; font-weight:500; border-bottom:1px solid #e2e8f0;">{escaped}</td>'
    elif any(w in col_lower for w in ["recover", "read"]):
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; white-space:nowrap;"><span style="display:inline-block; font-weight:600; color:#15803d; background:#ecfdf5; padding:3px 9px; border-radius:4px; border:1px solid #a7f3d0; font-size:12.5px;">{escaped}</span></td>'
    elif any(w in col_lower for w in ["action", "solution", "recommendation", "sheet"]):
        return f'<td style="padding:9px 12px; color:#0f766e; background:#f0fdfa; border-bottom:1px solid #e2e8f0;"><code style="color:#0f766e; background:none;">{escaped}</code></td>'
    else:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0;">{escaped}</td>'


def _build_html_table(rows: list, headers: list) -> str:
    th_html = "".join(
        f'<th style="padding:10px 12px; font-weight:600; color:#334155; border-bottom:2px solid #cbd5e1;">{html.escape(h)}</th>'
        for h in headers
    )
    tbody_rows = []
    for r in rows:
        tds = []
        for i, val in enumerate(r):
            h_name = headers[i] if i < len(headers) else ""
            tds.append(_format_cell(val, col_name=h_name, is_first=(i == 0)))
        tbody_rows.append(f'<tr>{"".join(tds)}</tr>')
    tbody_html = "".join(tbody_rows)
    return f"""<div style="overflow-x:auto; margin-top:10px; margin-bottom:12px;">
  <table class="table table-bordered table-striped" style="width:100%; border-collapse:collapse; background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; font-size:13.5px; text-align:left;">
    <thead><tr style="background:#f8fafc;">{th_html}</tr></thead>
    <tbody>{tbody_html}</tbody>
  </table>
</div>"""


def _parse_llm_json_to_table(parsed) -> tuple[list, list] | None:
    import re
    mapping_regex = re.compile(
        r"^([^:]+):\s*([A-Za-z0-9_\-+]+)\s*->\s*([A-Za-z0-9_\-+]+)(?:\s*\(([^)]+)\))?"
    )

    if isinstance(parsed, dict):
        for wrapper in ["problematic_lanes", "lanes", "results", "issues", "diagnoses", "items", "recoveries"]:
            if wrapper in parsed and isinstance(parsed[wrapper], (list, dict)):
                parsed = parsed[wrapper]
                break

    rows = []
    if isinstance(parsed, dict):
        for lane, details in parsed.items():
            if isinstance(details, dict):
                issue = ""
                action = ""
                other_details = []
                for k, v in details.items():
                    k_lower = k.lower()
                    if any(w in k_lower for w in ["issue", "problem", "hypothesis", "diagnosis"]):
                        issue = str(v)
                    elif any(w in k_lower for w in ["action", "solution", "recommendation", "sheet"]):
                        action = str(v)
                    else:
                        other_details.append(f"{k}: {v}")
                if not issue and details:
                    first_k, first_v = next(iter(details.items()))
                    issue = f"{first_k}: {first_v}"
                if other_details:
                    action = f"{action} ({'; '.join(other_details)})" if action else "; ".join(other_details)
                rows.append([str(lane), issue or "-", action or "-"])
            elif isinstance(details, (str, int, float)):
                rows.append([str(lane), str(details), "-"])
        if rows:
            return ["Lane / Item", "Identified Issue", "Recommended Action"], rows

    elif isinstance(parsed, list):
        # 1. Check if items represent unpacked/exploded mappings
        has_mappings = any(
            isinstance(item, dict) and any(
                "mapping" in k.lower() and isinstance(v, list) and len(v) > 0
                for k, v in item.items()
            )
            for item in parsed
        )

        if has_mappings:
            unpacked_rows = []
            for item in parsed:
                if not isinstance(item, dict):
                    continue
                issue = item.get("Issue", "-")
                rel = item.get("Relationship", issue)
                action = item.get("Action", item.get("RecommendedAction", "-"))
                gen_reads = item.get("RecoverableReads", "-")

                mappings = []
                for k, v in item.items():
                    if "mapping" in k.lower() and isinstance(v, list):
                        mappings = v
                        break

                if mappings:
                    for m_str in mappings:
                        match = mapping_regex.match(str(m_str).strip())
                        if match:
                            sid, exp_idx, obs_bc, rds = match.groups()
                            rds_clean = rds.replace("reads", "").replace("read", "").strip() if rds else ""
                            rds_str = f"{int(rds_clean):,}" if rds_clean.isdigit() else (rds or str(gen_reads))
                            unpacked_rows.append([sid, exp_idx, obs_bc, rel, rds_str, action])
                        else:
                            unpacked_rows.append([str(m_str), "-", "-", rel, str(gen_reads), action])
                else:
                    samples_list = item.get("AffectedSamples", item.get("Samples", []))
                    if isinstance(samples_list, list) and samples_list:
                        obs_val = ", ".join(str(x) for x in samples_list)
                    elif isinstance(samples_list, str) and samples_list:
                        obs_val = samples_list
                    else:
                        bc_match = re.search(r"[ACGTN]{6,}\+[ACGTN]{6,}", issue + " " + action)
                        obs_val = bc_match.group(0) if bc_match else "-"
                    unpacked_rows.append(["-", "-", obs_val, rel, str(gen_reads), action])

            if unpacked_rows:
                headers = [
                    "Candidate Sample",
                    "Expected Index",
                    "Found in Unassigned",
                    "Relationship / Similarity",
                    "Recoverable Reads",
                    "Recommended Action",
                ]
                return headers, unpacked_rows

        # 2. Check if list items represent direct mapping objects
        is_direct_mapping = any(
            isinstance(item, dict) and any(
                any(w in k.lower() for w in ["expectedindex", "observedbarcode", "unknownbarcode", "candidatesample"])
                for k in item.keys()
            )
            for item in parsed
        )
        if is_direct_mapping:
            unpacked_rows = []
            for item in parsed:
                if isinstance(item, dict):
                    sid = item.get("Sample", item.get("CandidateSample", item.get("SampleId", "-")))
                    exp_idx = item.get("ExpectedIndex", "-")
                    obs_bc = item.get("ObservedBarcode", item.get("UnknownBarcode", item.get("FoundInUnassigned", "-")))
                    rel = item.get("Relationship", item.get("Issue", "-"))
                    rds = item.get("RecoverableReads", item.get("Reads", "-"))
                    action = item.get("Action", item.get("RecommendedAction", "-"))
                    rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
                    unpacked_rows.append([str(sid), str(exp_idx), str(obs_bc), str(rel), rds_str, str(action)])
            if unpacked_rows:
                headers = [
                    "Candidate Sample",
                    "Expected Index",
                    "Found in Unassigned",
                    "Relationship / Similarity",
                    "Recoverable Reads",
                    "Recommended Action",
                ]
                return headers, unpacked_rows

        # 3. Check if list items represent grouped issues with sample lists
        is_grouped = any(
            isinstance(item, dict) and any(
                isinstance(v, list) for k, v in item.items() if any(w in k.lower() for w in ["sample", "barcode", "lane"])
            )
            for item in parsed
        )

        has_recovery = any(
            isinstance(item, dict) and any(
                any(w in k.lower() for w in ["recover", "read_count", "potential_reads", "reads"])
                for k in item.keys()
            )
            for item in parsed
        )

        if is_grouped:
            for item in parsed:
                if isinstance(item, dict):
                    issue = ""
                    action = ""
                    recovery = ""
                    samples = []
                    for k, v in item.items():
                        k_lower = k.lower()
                        if any(w in k_lower for w in ["issue", "problem", "hypothesis", "diagnosis"]):
                            issue = str(v)
                        elif any(w in k_lower for w in ["action", "solution", "recommendation", "sheet"]):
                            action = str(v)
                        elif any(w in k_lower for w in ["recover", "read_count", "potential_reads", "reads"]):
                            recovery = str(v)
                        elif any(w in k_lower for w in ["sample", "barcode", "lane"]):
                            if isinstance(v, list):
                                samples = v
                            elif isinstance(v, str):
                                samples = [v]
                    if has_recovery:
                        rows.append([issue or "-", recovery or "-", action or "-", samples if samples else "-"])
                    else:
                        rows.append([issue or "-", action or "-", samples if samples else "-"])
            if rows:
                headers = (
                    ["Identified Issue", "Potential Recoverable Reads", "Recommended Action", "Affected Samples / Barcodes"]
                    if has_recovery
                    else ["Identified Issue", "Recommended Action", "Affected Samples"]
                )
                return headers, rows
        else:
            for item in parsed:
                if isinstance(item, dict):
                    lane = ""
                    issue = ""
                    action = ""
                    recovery = ""
                    for k, v in item.items():
                        k_lower = k.lower()
                        if any(w in k_lower for w in ["lane", "sample", "id"]):
                            lane = str(v) if not isinstance(v, list) else ", ".join(str(x) for x in v)
                        elif any(w in k_lower for w in ["issue", "problem", "hypothesis", "diagnosis"]):
                            issue = str(v)
                        elif any(w in k_lower for w in ["action", "solution", "recommendation", "sheet"]):
                            action = str(v)
                        elif any(w in k_lower for w in ["recover", "read_count", "potential_reads", "reads"]):
                            recovery = str(v)
                    if not lane and item:
                        lane = f"Lane {len(rows)+1}"
                    if has_recovery:
                        rows.append([lane or "-", issue or "-", recovery or "-", action or "-"])
                    else:
                        rows.append([lane or "-", issue or "-", action or "-"])
                elif isinstance(item, str):
                    rows.append([f"Item {len(rows)+1}", str(item), "-"])
            if rows:
                headers = (
                    ["Lane / Item", "Identified Issue", "Potential Recoverable Reads", "Recommended Action"]
                    if has_recovery
                    else ["Lane / Item", "Identified Issue", "Recommended Action"]
                )
                return headers, rows

    return None


def generate_multiqc_html(
    content: str,
    output_path: str,
    sample_name: str = "Demultiplexing",
    input_data: dict = None,
):
    """Generate a MultiQC custom content HTML file."""
    clean_text = content.strip()
    if clean_text.startswith("```"):
        lines = clean_text.splitlines()
        if lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].startswith("```"):
            lines = lines[:-1]
        clean_text = "\n".join(lines).strip()

    # Parse LLM JSON if possible
    parsed = None
    try:
        parsed = json.loads(clean_text)
    except Exception:
        pass

    # Check if structured data is provided via input_data
    has_structured_data = bool(
        input_data
        and (
            input_data.get("CandidateRecoveries")
            or input_data.get("UnmatchedZeroReadSamples")
            or input_data.get("UnassignedAbundantBarcodes")
        )
    )

    if has_structured_data:
        headers = [
            "Candidate Sample",
            "Expected Index",
            "Found in Unassigned",
            "Relationship / Similarity",
            "Recoverable Reads",
            "Recommended Action",
        ]
        rows = []

        llm_unmatched_action = ""
        llm_unassigned_action = ""
        if parsed and isinstance(parsed, list):
            for item in parsed:
                if isinstance(item, dict):
                    iss = str(item.get("Issue", "")).lower()
                    act = str(item.get("Action", item.get("RecommendedAction", "")))
                    if any(w in iss for w in ["zero", "failure", "massive", "mismatch"]):
                        llm_unmatched_action = act
                    elif any(w in iss for w in ["unassigned", "unknown", "spike-in", "omitted"]):
                        llm_unassigned_action = act

        # 1. Candidate Recoveries
        for c in input_data.get("CandidateRecoveries", []):
            sid = c.get("CandidateSample", "-")
            lane = c.get("Lane", "")
            s_label = f"{sid} ({lane})" if lane else sid
            exp_idx = c.get("ExpectedIndex", "-")
            obs_bc = c.get("UnknownBarcode", "-")
            rel = c.get("Relationship", "-")
            rds = c.get("Reads", 0)
            rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
            action = c.get("RecommendedAction", "-")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        # 2. Unmatched Zero-Read Samples
        for s in input_data.get("UnmatchedZeroReadSamples", []):
            sid = s.get("SampleId", "-")
            lane = s.get("Lane", "")
            s_label = f"{sid} ({lane})" if lane else sid
            exp_idx = s.get("ExpectedIndex", "-")
            obs_bc = s.get("FoundInUnassigned", "None found")
            rel = s.get("Relationship", "Zero reads (no match in undetermined)")
            rds = s.get("Reads", 0)
            rds_str = str(rds)
            action = llm_unmatched_action or s.get("RecommendedAction", "Verify sample sheet index against library preparation design")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        # 3. Unassigned Abundant Barcodes
        for u in input_data.get("UnassignedAbundantBarcodes", []):
            lane = u.get("Lane", "")
            s_label = f"[Unassigned Barcode] ({lane})" if lane else "[Unassigned Barcode]"
            exp_idx = "-"
            obs_bc = u.get("UnknownBarcode", "-")
            rel = u.get("Relationship", "Abundant unassigned barcode")
            rds = u.get("Reads", 0)
            rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
            action = llm_unassigned_action or u.get("RecommendedAction", "Check for omitted sample or spike-in in library")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        is_ok = len(rows) == 0
        badge_color = "#28a745" if is_ok else "#dc3545"
        status_label = "PASS" if is_ok else "ATTENTION NEEDED"

        if is_ok:
            body_display = """<div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #22c55e; border-radius:4px; padding:12px 16px; margin:10px 0; color:#15803d; font-size:13.5px;">
  <strong>✓ Demultiplexing Evaluation:</strong> No index discrepancies or barcode assignment issues detected.
</div>"""
        else:
            summary_boxes_html = ""
            if parsed and isinstance(parsed, list):
                boxes = []
                for item in parsed:
                    if isinstance(item, dict):
                        iss = item.get("Issue", "")
                        act = item.get("Action", item.get("RecommendedAction", ""))
                        if iss or act:
                            boxes.append(f"""<div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:4px; padding:10px 14px; margin:8px 0; font-size:13px; color:#991b1b;">
  <strong>Identified Issue:</strong> {html.escape(str(iss))}<br>
  <strong>Recommendation:</strong> {html.escape(str(act))}
</div>""")
                if boxes:
                    summary_boxes_html = "\n".join(boxes)

            table_html = _build_html_table(rows, headers)
            body_display = f"{summary_boxes_html}\n{table_html}" if summary_boxes_html else table_html

    else:
        # Fallback when structured data is not available
        is_ok = (
            "no problems" in clean_text.lower()
            or "no issues" in clean_text.lower()
            or "all correct" in clean_text.lower()
        )
        badge_color = (
            "#28a745"
            if is_ok
            else "#e0a800"
            if "warn" in clean_text.lower()
            else "#dc3545"
        )
        status_label = "PASS" if is_ok else "ATTENTION NEEDED"
        escaped_content = html.escape(clean_text)

        if is_ok:
            body_display = """<div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #22c55e; border-radius:4px; padding:12px 16px; margin:10px 0; color:#15803d; font-size:13.5px;">
  <strong>✓ Demultiplexing Evaluation:</strong> No index discrepancies or barcode assignment issues detected.
</div>"""
        elif parsed is not None:
            table_data = _parse_llm_json_to_table(parsed)
            pretty_json = html.escape(json.dumps(parsed, indent=2))
            if table_data:
                headers, rows = table_data
                body_display = _build_html_table(rows, headers)
            else:
                body_display = f'<pre style="background:#f8f9fa; border:1px solid #e9ecef; border-radius:4px; padding:12px; margin-top:8px; font-size:13px; font-family:Menlo,Monaco,Consolas,monospace;"><code>{pretty_json}</code></pre>'
        else:
            body_display = f'<pre style="white-space:pre-wrap; background:#f8f9fa; border:1px solid #e9ecef; border-radius:4px; padding:12px; margin-top:8px; font-size:13px; font-family:Menlo,Monaco,Consolas,monospace;">{escaped_content}</pre>'

    html_content = f"""<!--
id: 'demultiplexing_llm_evaluation'
section_name: 'Demultiplexing LLM evaluation'
description: 'Automated LLM diagnosis of index assignments, undetermined reads, and demultiplexing statistics.'
-->
<div style="font-family: inherit; margin: 15px 0;">
  <div style="display: flex; align-items: center; margin-bottom: 8px;">
    <h4 style="margin: 0; font-size: 1.05em; font-weight: 600;">{html.escape(sample_name)}</h4>
    <span style="margin-left: 12px; padding: 2px 8px; border-radius: 4px; font-size: 0.78em; font-weight: bold; color: #fff; background-color: {badge_color};">
      {status_label}
    </span>
  </div>
  {body_display}
</div>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)


if __name__ == "__main__":
    main()
