import re

def replace_slashes_in_file(file_path, reverse=False):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    string_pattern = re.compile(r'(["\'])(.*?)(["\'])')
    cp_pattern = re.compile(r'\bcp\b ')

    with open(file_path, 'w') as file:
        for line in lines:
            def replace_slashes(match):
                string_content = match.group(2)
                if '\\r' in string_content or 'np.round' in string_content:
                    return match.group(0)
                if reverse:
                    string_content = string_content.replace('\\', '/')
                else:
                    string_content = string_content.replace('/', '\\\\')
                return f'{match.group(1)}{string_content}{match.group(3)}'

            new_line = string_pattern.sub(replace_slashes, line)
            new_line = cp_pattern.sub('copy ', new_line)
            file.write(new_line)

def process_files(file_paths, reverse=False):
    for file_path in file_paths:
        replace_slashes_in_file(file_path, reverse)

if __name__ == "__main__":
    files_to_process = [
        './Input_files/input_file_co2_only.txt',
        './Input_files/input_file_full_range.txt',
        './Input_files/input_file_short_range.txt',
        './utils.py',
        './retrieval-input.py',
        './plot_retrieval_results.py',
        './plot_mol_conditions-input.py',
        './plot_mol_contributions-input.py',
    ]
    reverse_option = False  # Set to True if you want to reverse the slashes
    process_files(files_to_process, reverse_option)
