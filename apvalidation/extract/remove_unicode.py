from unidecode import unidecode

def remove_unicode_from_paramter_file(file_path):
    print(f"removing unicode from {file_path}")
    # Open the file in binary mode
    with open(file_path, 'rb') as file:
        # Read the binary content
        binary_content = file.read()

    # Decode the binary content using a suitable encoding ('latin1')
    decoded_content = binary_content.decode('latin1', errors='replace')

    # Use unidecode to convert the content
    converted_content = unidecode(decoded_content)

    # Reopen the file in binary mode for writing (overwrite)
    with open(file_path, 'wb') as file:
        # Write the converted content back to the file
        file.write(converted_content.encode('utf-8'))