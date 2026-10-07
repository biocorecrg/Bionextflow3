#!/usr/bin/env python3
"""
prepare_demux_stats_for_llm.py

Filters and analyzes demultiplexing statistics (e.g. Illumina Stats.json or Aviti RunStats.json)
for consumption by LLM diagnosis services.

Key functions:
1. Extracts expected samples, reads, and top undetermined barcodes per lane.
2. Performs deterministic DNA sequence analysis between top undetermined barcodes
   and expected indices (testing i5/i7 reverse complement, swapped indices,
   Hamming distance/mismatches, and shifts).
3. Produces a structured JSON with explicit candidate recoveries and unassigned abundant barcodes.
"""

import argparse
import json
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter and prepare demux stats with index comparison for LLM diagnosis."
    )
    parser.add_argument(
        "--stats-json",
        "-i",
        required=True,
        help="Path to input Stats.json (Illumina) or RunStats.json (Aviti).",
    )
    parser.add_argument(
        "--output-json",
        "-o",
        default="demux_stats_llm.json",
        help="Path to output prepared JSON file (default: demux_stats_llm.json).",
    )
    parser.add_argument(
        "--top-barcodes",
        "-t",
        type=int,
        default=15,
        help="Number of top unknown barcodes to extract per lane (default: 15).",
    )
    parser.add_argument(
        "--max-mismatches",
        "-m",
        type=int,
        default=2,
        help="Maximum mismatches allowed when testing barcode similarities (default: 2).",
    )
    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    tr = str.maketrans("ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKWSVHBDNtgcaayrmkwsvhbdn")
    return seq.translate(tr)[::-1]


def split_index(index_seq: str) -> tuple[str, str]:
    """Split dual-index string into (i7, i5). For single-index, i5 is empty."""
    if not index_seq:
        return "", ""
    clean = index_seq.replace("-", "+")
    if "+" in clean:
        parts = clean.split("+", 1)
        return parts[0].strip(), parts[1].strip()
    return clean.strip(), ""


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings. Penalize length differences."""
    if not s1 and not s2:
        return 0
    if not s1 or not s2:
        return max(len(s1), len(s2))
    min_len = min(len(s1), len(s2))
    diff = sum(c1 != c2 for c1, c2 in zip(s1[:min_len], s2[:min_len]))
    diff += abs(len(s1) - len(s2))
    return diff


def compare_barcodes(obs_seq: str, exp_seq: str, max_mm: int = 2) -> dict | None:
    """
    Compare an observed unassigned barcode against an expected index.
    Tests hypotheses:
    - exact_match
    - i5_revcomp
    - i7_revcomp
    - swap (i7/i5 swapped)
    - swap_i5_revcomp
    - swap_both_revcomp
    - direct mismatch (Hamming distance)
    - i5_revcomp with mismatch
    """
    obs_i7, obs_i5 = split_index(obs_seq)
    exp_i7, exp_i5 = split_index(exp_seq)

    if not obs_i7 or not exp_i7:
        return None

    # Dual-index comparison
    if obs_i5 and exp_i5:
        exp_i5_rc = reverse_complement(exp_i5)
        exp_i7_rc = reverse_complement(exp_i7)

        # 1. Exact match (unexpectedly assigned to undetermined)
        if obs_i7 == exp_i7 and obs_i5 == exp_i5:
            return {
                "relationship": "Exact match (assigned to undetermined)",
                "mismatches": 0,
                "hypothesis": "exact_match",
                "action": "Check demultiplexing thresholds or sample sheet spelling",
            }

        # 2. i5 reverse complement
        d_i7 = hamming_distance(obs_i7, exp_i7)
        d_i5_rc = hamming_distance(obs_i5, exp_i5_rc)
        tot_rc_i5 = d_i7 + d_i5_rc
        if tot_rc_i5 <= max_mm:
            return {
                "relationship": f"i5 Reverse Complement ({tot_rc_i5} mismatches)" if tot_rc_i5 > 0 else "i5 Reverse Complement (Exact)",
                "mismatches": tot_rc_i5,
                "hypothesis": "i5_reverse_complement",
                "action": f"Update sample sheet i5 index to {obs_i5}",
            }

        # 3. i7 reverse complement
        d_i7_rc = hamming_distance(obs_i7, exp_i7_rc)
        d_i5 = hamming_distance(obs_i5, exp_i5)
        tot_rc_i7 = d_i7_rc + d_i5
        if tot_rc_i7 <= max_mm:
            return {
                "relationship": f"i7 Reverse Complement ({tot_rc_i7} mismatches)" if tot_rc_i7 > 0 else "i7 Reverse Complement (Exact)",
                "mismatches": tot_rc_i7,
                "hypothesis": "i7_reverse_complement",
                "action": f"Update sample sheet i7 index to {obs_i7}",
            }

        # 4. Swapped i7 and i5
        d_swp_i7 = hamming_distance(obs_i7, exp_i5)
        d_swp_i5 = hamming_distance(obs_i5, exp_i7)
        tot_swp = d_swp_i7 + d_swp_i5
        if tot_swp <= max_mm:
            return {
                "relationship": f"i7 and i5 Swapped ({tot_swp} mismatches)" if tot_swp > 0 else "i7 and i5 Swapped (Exact)",
                "mismatches": tot_swp,
                "hypothesis": "swap",
                "action": f"Swap i7 and i5 in sample sheet (i7={obs_i7}, i5={obs_i5})",
            }

        # 5. Swapped + i5 revcomp
        d_swp_rc = hamming_distance(obs_i7, exp_i5_rc) + hamming_distance(obs_i5, exp_i7)
        if d_swp_rc <= max_mm:
            return {
                "relationship": f"Swapped + i5 RevComp ({d_swp_rc} mismatches)",
                "mismatches": d_swp_rc,
                "hypothesis": "swap_and_rc",
                "action": f"Swap and reverse complement in sample sheet (i7={obs_i7}, i5={obs_i5})",
            }

        # 6. Direct mismatch (typo or sequencing error in standard orientation)
        tot_direct = d_i7 + d_i5
        if tot_direct <= max_mm:
            return {
                "relationship": f"{tot_direct}-base mismatch / typo",
                "mismatches": tot_direct,
                "hypothesis": "barcode_mismatch",
                "action": f"Correct index typo in sample sheet to {obs_seq} or allow --barcode-mismatches {tot_direct}",
            }

    # Single-index comparison
    else:
        exp_rc = reverse_complement(exp_i7)
        d_exact = hamming_distance(obs_i7, exp_i7)
        if d_exact == 0:
            return {"relationship": "Exact match", "mismatches": 0, "hypothesis": "exact_match", "action": "Check demux"}
        if d_exact <= max_mm:
            return {"relationship": f"{d_exact}-base mismatch", "mismatches": d_exact, "hypothesis": "barcode_mismatch", "action": f"Update index to {obs_i7}"}
        d_rc = hamming_distance(obs_i7, exp_rc)
        if d_rc <= max_mm:
            return {"relationship": f"Reverse Complement ({d_rc} mismatches)", "mismatches": d_rc, "hypothesis": "reverse_complement", "action": f"Reverse complement index to {obs_i7}"}

    return None


def process_illumina_stats(data: dict, top_n: int = 15, max_mm: int = 2) -> dict:
    """Process standard bcl2fastq/bcl-convert Stats.json."""
    expected_samples_by_lane = []
    unknown_barcodes_by_lane = []
    candidate_recoveries = []
    unassigned_abundant = []

    # Map lane -> list of samples
    conversion_results = data.get("ConversionResults", [])
    unknown_barcodes_raw = data.get("UnknownBarcodes", [])

    # Index unknown barcodes by lane
    unknown_map = {}
    for item in unknown_barcodes_raw:
        lane = item.get("Lane")
        barcodes = item.get("Barcodes", {})
        sorted_bcs = sorted(barcodes.items(), key=lambda x: x[1], reverse=True)[:top_n]
        unknown_map[lane] = sorted_bcs
        unknown_barcodes_by_lane.append({
            "Lane": lane,
            "TopBarcodes": dict(sorted_bcs),
        })

    for cr in conversion_results:
        lane = cr.get("LaneNumber")
        demux_results = cr.get("DemuxResults", [])
        lane_samples = []
        zero_or_low_samples = []

        for dr in demux_results:
            sample_id = dr.get("SampleId", "")
            reads = dr.get("NumberReads", 0)
            idx_metrics = dr.get("IndexMetrics", [])
            idx_seq = idx_metrics[0].get("IndexSequence", "") if idx_metrics else ""

            s_info = {
                "SampleId": sample_id,
                "Index": idx_seq,
                "Reads": reads,
            }
            lane_samples.append(s_info)
            if reads == 0:
                zero_or_low_samples.append(s_info)

        expected_samples_by_lane.append({
            "Lane": lane,
            "TotalSamples": len(lane_samples),
            "ZeroReadSamples": len(zero_or_low_samples),
        })

        # Compare top unknown barcodes against zero-read samples in this lane
        lane_unknowns = unknown_map.get(lane, [])
        for obs_bc, count in lane_unknowns:
            best_match = None
            best_sample = None

            for exp_sample in zero_or_low_samples:
                res = compare_barcodes(obs_bc, exp_sample["Index"], max_mm=max_mm)
                if res:
                    if best_match is None or res["mismatches"] < best_match["mismatches"]:
                        best_match = res
                        best_sample = exp_sample

            if best_match and best_sample:
                candidate_recoveries.append({
                    "Lane": f"Lane {lane}",
                    "UnknownBarcode": obs_bc,
                    "Reads": count,
                    "CandidateSample": best_sample["SampleId"],
                    "ExpectedIndex": best_sample["Index"],
                    "Relationship": best_match["relationship"],
                    "Mismatches": best_match["mismatches"],
                    "Hypothesis": best_match["hypothesis"],
                    "RecommendedAction": best_match["action"],
                })
            else:
                # If no match and high count (> 1M or top 3), flag as unassigned abundant
                if count >= 1_000_000 or (lane_unknowns and obs_bc == lane_unknowns[0][0]):
                    unassigned_abundant.append({
                        "Lane": lane,
                        "UnknownBarcode": obs_bc,
                        "Reads": count,
                        "Note": "No match found among expected samples (potential omitted sample or spike-in)",
                    })

    # Track unmatched zero read samples (samples with zero reads not in candidate_recoveries)
    recovered_samples_set = {
        (c.get("Lane"), c.get("CandidateSample")) for c in candidate_recoveries
    }

    unmatched_zero_map = {}
    for cr in conversion_results:
        lane = cr.get("LaneNumber")
        for dr in cr.get("DemuxResults", []):
            reads = dr.get("NumberReads", 0)
            if reads == 0:
                sid = dr.get("SampleId", "")
                if (f"Lane {lane}", sid) not in recovered_samples_set:
                    idx_m = dr.get("IndexMetrics", [])
                    idx_seq = idx_m[0].get("IndexSequence", "") if idx_m else ""
                    k = (sid, idx_seq)
                    if k not in unmatched_zero_map:
                        unmatched_zero_map[k] = []
                    unmatched_zero_map[k].append(lane)

    unmatched_zero_samples = []
    for (sid, idx_seq), lanes in unmatched_zero_map.items():
        lanes_sorted = sorted(lanes)
        if len(lanes_sorted) > 1 and lanes_sorted == list(range(lanes_sorted[0], lanes_sorted[-1] + 1)):
            lane_str = f"Lanes {lanes_sorted[0]}-{lanes_sorted[-1]}"
        elif len(lanes_sorted) > 1:
            lane_str = f"Lanes {', '.join(str(l) for l in lanes_sorted)}"
        else:
            lane_str = f"Lane {lanes_sorted[0]}" if lanes_sorted else "Unknown"

        unmatched_zero_samples.append({
            "SampleId": sid,
            "Lane": lane_str,
            "ExpectedIndex": idx_seq,
            "Reads": 0,
            "FoundInUnassigned": "None found",
            "Relationship": "Zero reads (no match in undetermined)",
            "RecommendedAction": "Verify sample sheet index against library preparation design",
        })

    # Consolidate unassigned abundant barcodes across lanes
    unassigned_map = {}
    for item in unassigned_abundant:
        bc = item["UnknownBarcode"]
        lane = item["Lane"]
        reads = item["Reads"]
        if bc not in unassigned_map:
            unassigned_map[bc] = {"lanes": [], "total_reads": 0, "notes": item.get("Note", "")}
        unassigned_map[bc]["lanes"].append(lane)
        unassigned_map[bc]["total_reads"] += reads

    consolidated_unassigned = []
    for bc, data in sorted(unassigned_map.items(), key=lambda x: x[1]["total_reads"], reverse=True)[:15]:
        lanes_sorted = sorted(data["lanes"])
        if len(lanes_sorted) > 1 and lanes_sorted == list(range(lanes_sorted[0], lanes_sorted[-1] + 1)):
            lane_str = f"Lanes {lanes_sorted[0]}-{lanes_sorted[-1]}"
        elif len(lanes_sorted) > 1:
            lane_str = f"Lanes {', '.join(str(l) for l in lanes_sorted)}"
        else:
            lane_str = f"Lane {lanes_sorted[0]}" if lanes_sorted else "Unknown"

        consolidated_unassigned.append({
            "UnknownBarcode": bc,
            "Lane": lane_str,
            "Reads": data["total_reads"],
            "FoundInUnassigned": bc,
            "Relationship": "Abundant unassigned barcode",
            "RecommendedAction": "Check for omitted sample or spike-in in library",
            "Note": data["notes"],
        })

    return {
        "ExpectedSamples": expected_samples_by_lane,
        "UnknownBarcodes": unknown_barcodes_by_lane,
        "CandidateRecoveries": candidate_recoveries,
        "UnmatchedZeroReadSamples": unmatched_zero_samples,
        "UnassignedAbundantBarcodes": consolidated_unassigned,
    }


def main():
    args = parse_args()

    try:
        with open(args.stats_json, "r", encoding="utf-8") as f:
            raw_data = json.load(f)
    except Exception as e:
        print(f"Error reading '{args.stats_json}': {e}", file=sys.stderr)
        sys.exit(1)

    # Process stats
    prepared = process_illumina_stats(
        raw_data,
        top_n=args.top_barcodes,
        max_mm=args.max_mismatches,
    )

    # Write output JSON
    try:
        with open(args.output_json, "w", encoding="utf-8") as f:
            json.dump(prepared, f, indent=2)
        print(f"Prepared demux stats saved to: {args.output_json}")
        print(f"  - Lanes analyzed: {len(prepared['ExpectedSamples'])}")
        print(f"  - Candidate recoveries identified: {len(prepared['CandidateRecoveries'])}")
        print(f"  - Unassigned abundant barcodes flagged: {len(prepared['UnassignedAbundantBarcodes'])}")
    except Exception as e:
        print(f"Error writing '{args.output_json}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
