import pandas as pd 
from path import Path
from jinja2 import Environment, FileSystemLoader  # html template engine
import visualize as bv



def generate_voc_html(feature, values, first_detected, world_map, state_map, genetic_distance_plot,
                  aa_distance_plot, s_aa_distance_plot):
    # express plots in html and JS
    world_map = plotly.offline.plot(world_map, include_plotlyjs=False, output_type='div')
    state_map = plotly.offline.plot(state_map, include_plotlyjs=False, output_type='div')
#     county_map = plotly.offline.plot(county_map, include_plotlyjs=False, output_type='div')
    genetic_distance_plot = plotly.offline.plot(genetic_distance_plot, include_plotlyjs=False, output_type='div')
    aa_distance_plot = plotly.offline.plot(aa_distance_plot, include_plotlyjs=False, output_type='div')
    s_aa_distance_plot = plotly.offline.plot(s_aa_distance_plot, include_plotlyjs=False, output_type='div')
    # generate output messages
    #TODO: expt_name, first_detected
    # dir containing our template
    file_loader = FileSystemLoader('templates')
    # load the environment
    env = Environment(loader=file_loader)
    # load the template
    template = env.get_template('voc.html')
    # render data in our template format
    html_output = template.render(feature=feature, values=values,
                                  world_map=world_map, state_map=state_map, 
                                  genetic_distance_plot=genetic_distance_plot, 
                                  aa_distance_plot=aa_distance_plot, s_aa_distance_plot=s_aa_distance_plot,
                                  first_detected=first_detected)
    return html_output


def generate_voc_data(feature, values, gisaid_data, tree_fp, subs_fp,
                  meta_fp, states_fp, patient_zero):
    genetic_distance_plot = bv.genetic_distance(tree_fp, meta_fp, patient_zero)
    aa_distance_plot = bv.aa_distance(subs_fp, meta_fp)
    s_aa_distance_plot = bv.s_aa_distance(subs_fp, meta_fp)
    state_map, _, _ = bv.map_by_state(gisaid_data, feature, values, states_fp)
    world_map, _, _ = bv.map_by_country(gisaid_data, feature, values)
    r = gisaid_data.loc[gisaid_data[feature].isin(values)]
    date = r['date_submitted'].min()
    state = r[r['date_submitted']==date]['division'].unique()
    cntry = r[r['date_submitted']==date]['country'].unique()
    first_detected = f"The {values} {feature} was first detected on {date} in {state}, {cntry}"
    return first_detected, genetic_distance_plot, aa_distance_plot, s_aa_distance_plot, world_map, state_map

def save_html(html_output: str, filename: str):
    with open(filename, 'w') as f:
        f.write(html_output)