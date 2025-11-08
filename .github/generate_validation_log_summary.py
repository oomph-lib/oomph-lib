import argparse
import os
import re


def print_cyan(text):
    (CYAN, RESET) = ("\033[96m", "\033[0m")
    print(f"{CYAN}{text}{RESET}")


def print_green(text):
    (GREEN, RESET) = ("\033[92m", "\033[0m")
    print(f"{GREEN}{text}{RESET}")


def print_orange(text):
    (ORANGE, RESET) = ("\033[93m", "\033[0m")
    print(f"{ORANGE}{text}{RESET}")


def print_red(text):
    (RED, RESET) = ("\033[91m", "\033[0m")
    print(f"{RED}{text}{RESET}")


def generate_markdown_summary(counts):
    summary = "### Validation Log Summary\n\n"
    summary += "| Status | Count |\n"
    summary += "|:-------|------:|\n"
    for (key, value) in counts.items():
        summary += f"| `{key}` | {value} |\n"
    summary += f"| **Total** | {sum(counts.values())} |\n"
    return summary


def count_occurrences(log_file, patterns):
    counts = {pattern: 0 for pattern in patterns}
    regex_patterns = re.compile("|".join(re.escape(p) for p in patterns))
    with open(log_file, "r", encoding="utf-8") as file:
        for line in file:
            matches = regex_patterns.findall(line)
            for match in matches:
                counts[match] += 1
    return counts


def parse_args():
    # fmt: off
    parser = argparse.ArgumentParser(description="Process a validation log file.")
    parser.add_argument("log_file", metavar="LOG_FILE", type=str, help="Path to the validation log file")
    parser.add_argument("--print-markdown-table", action="store_true", help="Path to the validation log file")
    args = parser.parse_args()
    # fmt: on
    if not os.path.exists(args.log_file):
        parser.error(f"Error: {args.log_file} not found.")
    return args


def main():
    args = parse_args()

    patterns = [
        ("[OK]", print_green),
        ("[FAILED]", print_red),
        ("[FAILED] I/O error", print_orange),
        ("Warning: files have different numbers of lines", print_orange),
    ]

    # Count the number of occurrences of each pattern listed above in the log file
    counts = count_occurrences(args.log_file, [p[0] for p in patterns])

    if args.print_markdown_table:
        # Print a Markdown table summary to standard output for redirection in GitHub Actions
        print(generate_markdown_summary(counts))
    else:
        # Just print a plain summary of the results
        print_cyan("\nTest results:")
        for (result, print_fn) in patterns:
            print_fn(f"  * {result:<50}: {counts[result]}")
        print()


if __name__ == "__main__":
    main()
