import sys

# Read all lines from stdin
names = sys.stdin.readlines()

# Strip newline characters and other whitespace from each name
names = [name.strip() for name in names]

# Join the names into a comma-separated string
comma_separated_names = ','.join(names)

# Output the result to stdout
sys.stdout.write(comma_separated_names + '\n')
