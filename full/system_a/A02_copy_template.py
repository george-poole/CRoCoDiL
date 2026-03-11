import os
import json

TEMPLATE_FILE_NAME = 'A01_template'

def copy_template(
    param_name: str,
    param_tex: str,
    file_name: str,
    heading: str,
    overwrite: bool,
    strip: bool,
) -> None:
    with open(f'{TEMPLATE_FILE_NAME}.ipynb', 'r') as f:
        nb = json.load(f)

    for cell in nb['cells']:
        if cell['cell_type'] == 'markdown':
            cell['source'] = [
                line if not line.startswith('# ') else f'# {heading}'
                for line in cell['source']
            ]
        if cell['cell_type'] == 'code':
            cell['source'] = [
                line if not line.startswith('PARAM_NAME = ') else f'PARAM_NAME = {repr(param_name)}'
                for line in cell['source']
            ]
            cell['source'] = [
                line if not line.startswith('PARAM_TEX = ') else f'PARAM_TEX = {repr(param_tex)}'
                for line in cell['source']
            ]
            if strip:
                cell['outputs'] = []
                cell['execution_count'] = None

    new_file_path = f'{file_name}.ipynb'
    if not overwrite and os.path.exists(new_file_path):
        print(f'Not overwriting {new_file_path}')
    else:
        with open(new_file_path, 'w') as f:
            json.dump(nb, f, indent=1)


if __name__ == '__main__':
    OVERWRITE = False
    STRIP = True
    param_name = ('Ra', 'Da', 'sr')
    param_tex = ('Ra', 'Da', 's_r')
    file_name = ('A11_rayleigh', 'A12_damkohler', 'A13_saturation')
    heading = ('Rayleigh numbers', 'Damköhler numbers', 'Saturations')
    for args in zip(
        param_name, param_tex, file_name, heading
    ):
        copy_template(*args, OVERWRITE, STRIP)
