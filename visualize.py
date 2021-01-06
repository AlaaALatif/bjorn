import pandas as pd
from path import Path
import plotly
import plotly.express as px
import plotly.graph_objects as go
import onion_trees as bv
from urllib.request import urlopen
import json
import statsmodels as sm
from statsmodels.formula.api import ols
import mutations as bm


def aa_distance(subs_fp, meta_fp, alpha=0.05):
    alab_subs = pd.read_csv(subs_fp)
    alab_subs['nonsyn'] = False
    alab_subs.loc[alab_subs['ref_aa']!=alab_subs['alt_aa'], 'nonsyn'] = True
    alab_subs['S_nonsyn'] = False
    alab_subs.loc[(alab_subs['gene']=='S') & (alab_subs['ref_aa']!=alab_subs['alt_aa']), 'S_nonsyn'] = True
    dists_df = (alab_subs.groupby('fasta_hdr')
                .agg(num_nonsyn_muts=('nonsyn', 'sum'), num_S_nonsyn_muts=('S_nonsyn', 'sum'))
                .reset_index())
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df['date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df['month'] = df['date'].dt.month
    df['doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()]
    model = ols('num_nonsyn_muts ~ doy', data=df).fit()
    df['predict'] = model.predict(df['doy'])
    df['p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df['outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(data=go.Scatter(y=df[df['outlier']==False]['num_nonsyn_muts'], x=df[df['outlier']==False]['doy'], 
                                name='samples', mode='markers', marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['num_nonsyn_muts'], x=df[df['outlier']==True]['doy'],
                             mode='markers', marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['doy'], name='OLS', mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Amino Acid Changes (root-to-tip)', xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def s_aa_distance(subs_fp, meta_fp, alpha=0.05):
    alab_subs = pd.read_csv(subs_fp)
    alab_subs['nonsyn'] = False
    alab_subs.loc[alab_subs['ref_aa']!=alab_subs['alt_aa'], 'nonsyn'] = True
    alab_subs['S_nonsyn'] = False
    alab_subs.loc[(alab_subs['gene']=='S') & (alab_subs['ref_aa']!=alab_subs['alt_aa']), 'S_nonsyn'] = True
    dists_df = (alab_subs.groupby('fasta_hdr')
                .agg(num_nonsyn_muts=('nonsyn', 'sum'), num_S_nonsyn_muts=('S_nonsyn', 'sum'))
                .reset_index())
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df['date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df['month'] = df['date'].dt.month
    df['doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()]
    model = ols('num_S_nonsyn_muts ~ doy', data=df).fit()
    df['predict'] = model.predict(df['doy'])
    df['p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df['outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(data=go.Scatter(y=df[df['outlier']==False]['num_S_nonsyn_muts'], x=df[df['outlier']==False]['doy'], 
                                name='samples', mode='markers', marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['num_S_nonsyn_muts'], x=df[df['outlier']==True]['doy'],
                             mode='markers', marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['doy'], name='OLS', mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Amino Acid Changes in the S protein(root-to-tip)', xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def genetic_distance(tree_fp, meta_fp, patient_zero, alpha=0.05):
    tree = bv.load_tree(tree_fp, patient_zero)
    dists = {n.name: tree.distance(n.name, patient_zero) for n in tree.get_terminals()}
    dists_df = (pd.DataFrame(index=dists.keys(), data=dists.values(), 
                      columns=['genetic_distance'])
         .reset_index()
         .rename(columns={'index': 'fasta_hdr'}))
    meta = pd.read_csv(meta_fp)
    sd_meta = meta[meta['location'].str.contains('San Diego')]
    df = pd.merge(dists_df, sd_meta, on='fasta_hdr')
    df['date'] = pd.to_datetime(df['collection_date'], errors='coerce')
    df['month'] = df['date'].dt.month
    df['doy'] = df['date'].dt.dayofyear
    df = df.loc[~df['doy'].isna()]
    model = ols('genetic_distance ~ doy', data=df).fit()
    df['predict'] = model.predict(df['doy'])
    df['p'] = model.outlier_test(method='fdr_bh')['fdr_bh(p)']
    df['outlier'] = False
    df.loc[df['p']<alpha, 'outlier'] = True
    fig = go.Figure(data=go.Scatter(y=df[df['outlier']==False]['genetic_distance'], x=df[df['outlier']==False]['doy'], 
                                name='samples', mode='markers', marker_color='rgba(30,144,255,.6)'))
    fig.add_trace(go.Scatter(y=df[df['outlier']==True]['genetic_distance'], x=df[df['outlier']==True]['doy'],
                             mode='markers', marker_color='rgba(220,20,60,.6)', name='SoIs',
                 text=df[df['outlier']==True][['ID', 'date']],
                 hovertemplate = 
                 "<b>%{text[0]}</b><br>" +
                 "<b>%{text[1]}</b><br>"))
    fig.add_trace(go.Scatter(y=df['predict'], x=df['doy'], name='OLS', mode='lines', line_color='rgba(0,0,0,1.)'))
    fig.update_layout(yaxis_title='Genetic Distance (root-to-tip)', xaxis_title='Collection Date',
                      template='plotly_white', autosize=True)#, height=850, width=800)
    return fig


def map_by_state(data: pd.DataFrame, feature: str, values: list, states_fp: str):
    with open(states_fp) as f:
        states = json.load(f)
    state_map = {x['properties']['name']: x['id'] for x in states['features']}
    results = data.loc[(data[feature].isin(values)) & (data['country']=='United States of America')]
    results_by_state = results.groupby('division').agg(num_samples=('idx', 'nunique')).reset_index()
    results_by_state['id'] = results_by_state['division'].apply(lambda x: state_map.get(x, 'unk'))
#     fig = px.choropleth(results_by_state, geojson=states, scope="usa",
#                                locations='id', color='num_samples',# locationmode='USA-states',
#                                color_continuous_scale="bluered",
#                                range_color=(0, results_by_state['num_samples'].max()),
# #                                labels={'num_samples': f'Number of samples with {values}: ', 'division': 'loc:'},
#                                hover_data=['division', 'num_samples']
#                               )
    fig = px.choropleth_mapbox(results_by_state, geojson=states, 
                               locations='id', color='num_samples',
                               color_continuous_scale="bluered", center={"lat": 37.0902, "lon": -95.7129},
                               range_color=(0, results_by_state['num_samples'].max()),
                               mapbox_style="carto-positron", zoom=3,
                               opacity=0.5,
                               hover_data=['division', 'num_samples']
                               #labels={'num_samples':f'Number of samples with {values}', 'division': f'location:'}
                              )
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    return fig, state_map, results_by_state


def map_by_country(data: pd.DataFrame, feature: str, values: list):
    with urlopen('https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json') as response:
        countries = json.load(response)
    for c in countries['features']:
        if c['id']=='USA':
            assert c['properties']['name'] == 'United States of America'
    country_map = {x['properties']['name']: x['id'] for x in countries['features']}
    results = data.loc[data[feature].isin(values)]
    results_by_cntry = results.groupby('country').agg(num_samples=('idx', 'nunique')).reset_index()
    results_by_cntry['id'] = results_by_cntry['country'].apply(lambda x: country_map.get(x, 'unk'))
    fig = px.choropleth_mapbox(results_by_cntry, geojson=countries, 
                               locations='id', color='num_samples',
                               color_continuous_scale="bluered",
                               range_color=(0, results_by_cntry['num_samples'].max()),
                               mapbox_style="carto-positron", zoom=1,
                               opacity=0.5,
                               labels={'num_samples':f'Number of samples with {values}', 'id': f'location:'}
                              )
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    return fig, country_map, results_by_cntry