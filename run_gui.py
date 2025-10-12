from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import numpy as np
from helper_functions import *

app = Dash(__name__)

# Layout
app.layout = html.Div([
    html.H1("3D Helix: adjust pitch and radius"),

    html.Div([
        html.Label("Outer Length of each segment "),
        dcc.Input(id='outer-length', type='number', value=30, min=0.1, max=1000, step=0.1),
        html.Span(" mm")
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px'}),

    html.Div([
        html.Label("cut angle "),
        dcc.Input(id='cut-angle', type='number', value=7.5, min=0.1, max=89.9, step=0.1),
        html.Span(" degrees")
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

    html.Div([
        html.Label("pipe diameter "),
        dcc.Input(id='pipe-diameter', type='number', value=2.5, min=0.1, max=20, step=0.1),
        html.Span(" inches")
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

    html.Div([
        html.Label("twist angle "),
        dcc.Input(id='twist-angle', type='number', value=15, min=0.1, max=179, step=0.1),
        #dcc.Slider(id = 'twist-angle', value = 15, min = 0.1, max = 89, step = 0.1),
        html.Span(" degrees")
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

    html.Div([
        html.Label("desired length of helix "),
        dcc.Input(id='length', type='number', value = 1000, min=0.1, max=2000, step=0.1),
        html.Span(" mm")
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

    html.Div([
        html.Label("Number of segments"),
        dcc.Input(id='segments', type='number', value=45, min=4, max=200, step=1),
    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

    html.Div([
        dcc.Checklist(
            id='double_helix',
            options=[{'label': 'Display double helix', 'value': 'True'}],
            value=[],  # empty = unchecked
            inline=True
    ),    ], style={'display': 'block', 'alignItems': 'center', 'gap': '15px', 'marginTop': '15px'}),

       # Right-side figure
    html.Div([
        dcc.Graph(id='helix-plot', style={'height': '100%', 'width': '100%'})
    ], style={
        'position': 'fixed',
        'right': '0',
        'top': '0',
        'width': '50%',
        'height': '100%',
        'aspectRatio': '2 / 1',  # CSS aspect ratio
        'display': 'flex',
        'alignItems': 'center',
        'justifyContent': 'center',
        'backgroundColor': 'white'
    }),



html.Div(id='radius', style={'marginTop': '10px'}),
html.Div(id='actual_length', style={'marginTop': '10px'}),
html.Div(id='deltaz', style={'marginTop': '10px'}),
html.Div(id='segments_from_length', style={'marginTop': '10px'}),
html.Div(id='actual_length_from_length', style={'marginTop': '10px'}),
html.Div(id='L', style={'marginTop': '10px'}),
html.Div(id='Nturns', style={'marginTop': '10px'}),
html.Div(id='Nturns_from_length', style={'marginTop': '10px'}),



    ],
    style = {'display': 'block'}
    ),







# --- Update helix plot ---
@app.callback(
    Output('helix-plot', 'figure'),
    Output('radius', 'children'),
    Output('actual_length', 'children'),
    Output('deltaz', 'children'),
    Output('segments_from_length', 'children'),
    Output('actual_length_from_length', 'children'),
    Output('L', 'children'),
    Output('Nturns', 'children'),
    Output('Nturns_from_length', 'children'),

    
    Input('outer-length', 'value'),
    Input('cut-angle', 'value'),
    Input('pipe-diameter', 'value'),
    Input('twist-angle','value'),
    Input('length', 'value'),
    Input('segments', 'value'),
    Input('double_helix', 'value'),
)


def update_helix(outer_length, cut_angle, pipe_diameter, twist_angle, length, segments, double_helix = ['True']):
    

    pipe_diameter *= 25.4 # mm to inches

    L = center_length_from_outer(outer_length, cut_angle, pipe_diameter)
    theta = cut_angle * 2
    tau = twist_angle

    if 'True' in double_helix:
        double_helix = True

    points_straight, helix_radius, _, actual_length, deltaz, segments_from_length, actual_length_from_length = full_helix_calculation(L, theta, tau, length = length, segments = segments)

    polar_angle = polar_angle_between(points_straight[1], points_straight[0])

    Nturns = find_number_of_turns(polar_angle, segments)
    Nturns_from_length = find_number_of_turns(polar_angle, segments_from_length)

    x, y, z = points_straight.T


    fig = go.Figure()

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, line=dict(width=20)))

    if double_helix:
        x, y, z = turn_spiral_about_z_axis(points_straight).T

        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, line=dict(width=20)))



    fig.update_layout(
        template='plotly_white',
        scene=dict(aspectmode='data'),
    )

    return fig,     \
    f"Helix diamter: {helix_radius * 2:.2f} mm",\
    f"Actual length using {segments} segments: {actual_length:.2f} mm",\
    f"height increase per segment: {deltaz:.2f} mm",\
    f"Number of segments to reach desired length ({length:.2f}): {segments_from_length}",\
    f"Actual length using {segments_from_length} segments: {actual_length_from_length:.2f} mm",\
    f"center length of the pipe: {L:.2f} mm",\
    f"Number of turns (calculated for {segments} segments): {Nturns:.2f}",\
    f"N turns (calculated for {segments_from_length} segments): {Nturns_from_length:.2f}"\

if __name__ == "__main__":
    app.run(debug=False)
