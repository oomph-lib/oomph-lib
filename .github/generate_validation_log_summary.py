import argparse
import os
import re


def count_occurrences(log_file, patterns):
    counts = {pattern: 0 for pattern in patterns}
    regex_patterns = re.compile("|".join(re.escape(p) for p in patterns))
    with open(log_file, "r", encoding="utf-8") as file:
        for line in file:
            matches = regex_patterns.findall(line)
            for match in matches:
                counts[match] += 1
    return counts


def generate_markdown_summary(counts):
    summary = "### Validation Log Summary\n\n"
    summary += "| Status | Count |\n"
    summary += "|:-------|------:|\n"
    for (key, value) in counts.items():
        summary += f"| `{key}` | {value} |\n"
    summary += f"| **Total** | {sum(counts.values())} |\n"
    return summary


def main():
    parser = argparse.ArgumentParser(description="Process a validation log file.")
    parser.add_argument("log_file", type=str, help="Path to the validation log file")
    args = parser.parse_args()

    if not os.path.exists(args.log_file):
        print(f"Error: {args.log_file} not found.")
        return

    patterns = [
        "[OK]",
        "[FAILED]",
        "[FAILED] I/O error",
        "Warning: files have different numbers of lines"
    ]

    counts = count_occurrences(args.log_file, patterns)
    summary = generate_markdown_summary(counts)

    # Print summary to standard output for redirection in GitHub Actions
    print(summary)


if __name__ == "__main__":
    main()
