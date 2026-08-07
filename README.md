# Double helix from pipe segments

A small project that calculates optimal parameters of a intertwined doulbe helix made of pipe segments cut at an angle.

See the web version at [bojanovsky.at/helix](https://bojanovsky.at/helix).


## Details
The testing environments of the calculation can be seen in the [Python Notebook](Joint_face_normal_method.ipynb).
All the finished and somewhat working functions are in the file [helper_funcitons.py](helper_functions.py).

To visualize and correctly compute the spiral, run run_gui.py in this directory
    
    python run_gui.py

To run this, you will need numpy, plotly and dash

    pip install numpy plotly dash

Clone this directory using

    git clone https://github.com/jocelas/spiral

The theory is mostly adapted from [https://github.com/PubInv/segmented-helixes/](https://github.com/PubInv/segmented-helixes/), but adapted to our case.

For things to do, look at the issues on the [github page](https://github.com/jocelas/spiral) of this project.

