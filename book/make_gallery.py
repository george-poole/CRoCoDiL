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
    captions: Iterable[str],
    definitions: Iterable[tuple[str, str]],
    overviews: Iterable[tuple[str, str]],
):
    markdown_string = heading(topheading, 1)

    for n, cpt, defn, ovvw in zip(subheadings, captions, definitions, overviews, strict=True):
        markdown_string = '\n'.join((markdown_string, heading(n, 2)))
        markdown_string = '\n'.join((markdown_string, cpt))
        td_schematic = table_data(*defn, html=html)
        td_overview = table_data(*ovvw, html=html)
        tds = (td_schematic, td_overview)
        Ntr = len(tds)
        tds_divided = [tds[i * Ntr: Ntr * (i + 1)] for i in range(ceil(len(tds) / Ntr))]
        trs = [table_row(tdsd) for tdsd in tds_divided if tdsd]
        markdown_string = '\n'.join((markdown_string, table(trs)))

    with open(gallery_file, "w") as f:
        f.write(markdown_string)


if __name__ == "__main__":
    degreek = lambda s: s.replace('θ', 'Theta').replace('μ', 'Mu')
    get_dir = lambda s: degreek(s.lower().replace(' ', '_'))
    get_letter = lambda s: degreek(s.split(' ')[1])
    HEADING = 'Gallery'
    SUBHEADINGS = (
        'System A', 
        'System Aμ',
        'System Aθ',  
        'System B', 
        'System C', 
        'System D',
    )
    CAPTIONS = (
        'Solutal convective dissolution in a rectangle',
        'Viscous fingering effects on convective dissolution in a rectangle',
        'Thermal buoyancy effects on convective dissolution in a rectangle',
        'Solutal convective dissolution in an anticline with a background flow',
        'Exchange flow dissolution in an inclined rectangle',
        'Buoyant plumes in thermosolutal convective dissolution',
    )
    SCHEMATICS = tuple(
        (
            f'./{get_dir(s)}/expo/{get_letter(s)}00_definition.html', 
            f'./{get_dir(s)}/expo/figures/{get_dir(s)}_Omega_0.png',
            'Definition',
        )
        for s in SUBHEADINGS
    )
    OVERVIEWS = tuple(
        (
            f'./{get_dir(s)}/expo/{get_letter(s)}10_overview.html', 
            f'./{get_dir(s)}/expo/figures/{get_letter(s)}10_overview/thumbnail.png',
            'Overview',
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
            CAPTIONS,
            SCHEMATICS, 
            OVERVIEWS, 
        )