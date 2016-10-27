
# Extract code blocks from an .Rmd file

import sys, textwrap

print "# This file is generated from the corresponding .Rmd file"
print
print

in_code = False
show_code = True
in_challenge = False
challenge_level = 0
for line in sys.stdin:
    line = line.rstrip()
    if line.startswith("```"):
        print
        in_code = not in_code
        if in_code:
            show_code = "echo=FALSE" not in line
        assert in_code or line == "```", line
    elif in_code:
        if show_code:
            print line
    elif line.startswith("#"):
        print "#" if in_challenge else ""
        n = line.count("#")
        if not in_challenge or n <= challenge_level:
            in_challenge = "{.challenge}" in line
            if in_challenge:
                challenge_level = n
        if in_challenge:
            line = line.replace("{.challenge}","").rstrip()
        banner = "#"*n + " " + (" " if n > 3 else ("-" if n > 2 else "=")) * (len(line)-n-1)
        print banner
        print line
        print banner
    elif in_challenge:
        for line2 in textwrap.wrap(line) or [""]:
            print "# " + line2
