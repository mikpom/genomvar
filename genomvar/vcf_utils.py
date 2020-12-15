import os
from jinja2 import Environment,FileSystemLoader

tmpl_dir = os.path.join(os.path.dirname(__file__),'tmpl')
env = Environment(
    loader=FileSystemLoader(tmpl_dir),
)
header_simple = env.get_template('vcf_head_simple.tmpl')
header = env.get_template('vcf_head.tmpl')
row_tmpl = env.get_template('vcf_row.tmpl')
