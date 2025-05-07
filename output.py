import re

# Function to mark start codons in an RTF file
def mark_start_codons(input_rtf, output_rtf):
    with open(input_rtf, 'r') as file:
        content = file.read()

    # Regular expression to find 'ATG' codons
    marked_content = re.sub(r'ATG', r'\b\{\\b **ATG**\\b\}', content)

    # Save the modified content to a new RTF file
    with open(output_rtf, 'w') as file:
        file.write(marked_content)

# Provide the input and output file paths
input_rtf = 'mitfa.rtf'
output_rtf = 'output_file_marked.rtf'

# Run the function
mark_start_codons(input_rtf, output_rtf)

