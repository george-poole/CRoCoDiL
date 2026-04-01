import json
from math import ceil
from typing import Iterable


def heading(
    text: str, 
    depth: int,
    ):

    return f"{'#' * depth} {text}\n"

def table(
    table_rows: list[str],
) -> str:
    if not table_rows:
        return ''
    br = "\n"
    return f"""<table>
    {br.join(i for i in table_rows)}
    </table>{br}"""


def table_row(
    table_datas: list[str],
) -> str:  
    if not table_datas:
        return ''
    br = "\n"
    return f"""<tr>
    {br.join(i for i in table_datas)}
    </tr>"""


def table_data(
    caption_path: str,
    fig_path: str,
    caption: str = '',
    img_width: int = 400,
    html: bool = False,
) -> str:
    if html: 
        caption_path = caption_path.replace('ipynb', 'html')
    return f"""<td>
    <img src="{fig_path}" width="{img_width}"/>
    <a href="{caption_path}">{caption}</a> 
    </td>"""


def extract_heading(
    file_path: str,
    default: str | None = None,
) -> str | None:
    
    with open(file_path, 'r', encoding='utf-8') as f:
        notebook: dict[str, list[dict]] = json.load(f)

    for cell in notebook.get('cells', []):
        if cell.get('cell_type') == 'markdown':
            for line in cell.get('source', []):
                stripped = line.strip()
                if stripped.startswith('# '):
                    return stripped[2:].strip()
                
    return default


def make_gallery(
    gallery_file: str,
    html: bool,
    topheading: str,
    subheadings: Iterable[str],
    schematics: Iterable[tuple[str, str]],
    overviews: Iterable[tuple[str, str]],
):
    markdown_string = heading(topheading, 1)

    for n, schm, ovvw in zip(subheadings, schematics, overviews, strict=True):
        markdown_string = '\n'.join((markdown_string, heading(n, 2)))
        td_schematic = table_data(*schm, html=html)
        td_overview = table_data(*ovvw, html=html)
        tds = (td_schematic, td_overview)
        Ntr = len(tds)
        tds_divided = [tds[i * Ntr: Ntr * (i + 1)] for i in range(ceil(len(tds) / Ntr))]
        trs = [table_row(tdsd) for tdsd in tds_divided if tdsd]
        markdown_string = '\n'.join((markdown_string, table(trs)))

    with open(gallery_file, "w") as f:
        f.write(markdown_string)


if __name__ == "__main__":
    HEADING = 'Gallery'
    SUBHEADINGS = (
        'System A', 'System B', 'System C', 'System D',
    )
    get_dir = lambda s: s.lower().replace(' ', '_')
    get_letter = lambda s: s.split(' ')[1]
    SCHEMATICS = tuple(
        (
            f'./{get_dir(s)}/expo/{get_letter(s)}00_definition.md', 
            f'./{get_dir(s)}/expo/figures/{get_dir(s)}_Omega_0.png',
        )
        for s in SUBHEADINGS
    )
    OVERVIEWS = tuple(
        (
            f'./{get_dir(s)}/expo/{get_letter(s)}10_overview.ipynb', 
            f'./{get_dir(s)}/expo/figures/{get_letter(s)}10_overview/thumbnail.png',
        )
        for s in SUBHEADINGS
    )
    OUTPUTS = (
        ('gallery.md', True), 
        ('gallery_local.md', False),
    )
    for file, html in OUTPUTS:
        make_gallery(
            file, 
            html,
            HEADING, 
            SUBHEADINGS, 
            SCHEMATICS, 
            OVERVIEWS, 
        )