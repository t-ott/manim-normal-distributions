# Animating Normal Distributions with ```manim```

Python code to render two and three-dimensional scenes of normal distributions with the ```manim``` library, specifically with the [Manim Community](https://github.com/ManimCommunity/manim) fork of ```manim```.

***
**A blog post walking through some of the code can be found [here](https://www.t-ott.dev/2021/11/24/animating-normal-distributions).**
***

### Univariate Normal Distributions
The file ```univariate.py``` has ```manim``` scenes which animate adjustments to the mean and standard deviation of a univariate normal distribution.

### Bivariate Normal Distributions
The file ```bivariate.py``` has ```manim``` scenes using ```ThreeDScene``` which animate adjustments to the means, standard deviations, and correlation factor of bivariate normal distributions.

***

### Rendering Scenes
If ```manim``` and its various dependencies are installed (see ```manim``` [docs](https://docs.manim.community/en/stable/installation.html) for installation instructions), you can render these scenes by navigating to the project directory and running a command such as:
```
manim -pqm bivariate.py AdjustRho
```
Which would render the ```AdjustRho``` scene. The ```-p``` flag would play the rendered video once it's complete, and ```-qm``` specifies a medium quality render (720p).
