#!/bin/bash

# Script to convert a text file (like bare_sample.txt) to a basic HTML page.
# Usage: ./txt_to_html.sh <input_text_file>
# Example: ./txt_to_html.sh /home/wilkie/code/pidp11-astrology/bare_sample.txt

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_text_file>" >&2
    exit 1
fi

input_file="$1"

if [ ! -r "$input_file" ]; then
    echo "Error: Input file '$input_file' not found or is not readable." >&2
    exit 1
fi

# HTML Header
cat <<'EOF'
<html>
<head>
<title>PiDP Astrological Wisdom</title>
</head>
<body>
EOF

# Process the input file
while IFS= read -r line; do
  case "$line" in
    "Calculating for date:"*)
      echo "<h1>$line</h1>"
      ;;
    "  Wisdom for "*) # Matches lines starting with two spaces then "Wisdom for "
      content=$(sed 's/^[[:blank:]]*//' <<< "$line")
      echo "<h2>$content</h2>"
      ;;
    "    "*) # Matches lines starting with four spaces
      content=$(sed 's/^[[:blank:]]*//' <<< "$line")
      echo "<p>$content</p>"
      ;;
  esac
done < "$input_file"

# HTML Footer
cat <<'EOF'
</body>
</html>
EOF
