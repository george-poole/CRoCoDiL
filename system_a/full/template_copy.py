import os
import json

def copy_template(
    template_file: str,
    param_name: str,
    param_tex: str,
    file_name: str,
    heading: str,
    overwrite: bool,
    outputs: bool,
) -> None:
    with open(f'{template_file}.ipynb', 'r') as f:
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
            if not outputs:
                cell['outputs'] = []
                cell['execution_count'] = None

    new_file_path = f'{file_name}.ipynb'
    if not overwrite and os.path.exists(new_file_path):
        print(f'Not overwriting {new_file_path}')
    else:
        with open(new_file_path, 'w') as f:
            json.dump(nb, f, indent=1)


if __name__ == '__main__':
    TEMPLATE_FILE_NAME = 'A01_template'
    OVERWRITE = True
    OUTPUTS = False
    EXECUTE = False
    param_name = ('Ra', 'Da', 'sr')
    param_tex = ('Ra', 'Da', 's_r')
    file_name = ('A11_rayleigh', 'A12_damkohler', 'A13_saturation')
    heading = ('Rayleigh numbers', 'Damköhler numbers', 'Residual saturations')
    for args in zip(
        param_name, param_tex, file_name, heading
    ):
        copy_template(TEMPLATE_FILE_NAME, *args, OVERWRITE, OUTPUTS)

    if EXECUTE:
        for fn in file_name:
            os.system(f"jupyter nbconvert --execute --to notebook --inplace {fn}.ipynb")
