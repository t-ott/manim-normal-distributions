from manim import *
import math


def PDF_normal(x, mu, sigma):
    '''
    General form of probability density function of univariate normal distribution
    '''
    return math.exp(-((x-mu)**2)/(2*sigma**2))/(sigma*math.sqrt(2*math.pi))


class AdjustMu(Scene):
    '''
    Scene to observe how adjustments to the mean of a normal distrubtion
    influences the shape of its probability density function
    '''

    def construct(self):
        ax = Axes(
            x_range = [-5, 5, 1],
            y_range = [0, 0.5, 0.1],
            axis_config = {'include_numbers':True}
        )

        # Initialize mu (distribution mean) ValueTracker to 0
        mu = ValueTracker(0)

        # Text to display distrubtion mean
        mu_text = MathTex(r'\mu = ').next_to(ax, UP, buff=0.2).set_color(YELLOW)
        # Always redraw the decimal value for mu for each frame
        mu_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(mu.get_value())
            .next_to(mu_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        # Define PDF curve, always redraw for each frame
        curve = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), 1), color=YELLOW)
        )

        # Start animation
        self.add(ax, mu_text, mu_value_text)
        self.play(Create(curve))
        self.play(
            mu.animate.set_value(2), run_time=1,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu.animate.set_value(-2), run_time=1.5,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu.animate.set_value(0), run_time=1,
            rate_func=rate_functions.smooth
        )
        self.play(Uncreate(curve))


class AdjustSigma(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distrubtion influences the shape of its probability density function
    '''

    def construct(self):
        ax = Axes(
            x_range = [-4, 4, 1],
            y_range = [0, 1, 0.2],
            axis_config = {'include_numbers':True}
        )

        # Initialize sigma (distribution standard deviation) ValueTracker to 1
        sigma = ValueTracker(1)

        # Text to display distrubtion standard deviation
        sigma_text = MathTex(r'\sigma = ').next_to(ax, UP, buff=0.2).set_color(YELLOW)
        # Always redraw the decimal value for sigma for each frame
        sigma_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(sigma.get_value())
            .next_to(sigma_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        # Define PDF curve, always redraw for each frame
        curve = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, 0, sigma.get_value()), color=YELLOW)
        )

        # Start animation
        self.add(ax, sigma_text, sigma_value_text)
        self.play(Create(curve))
        self.play(
            sigma.animate.set_value(1.5), run_time=1,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            sigma.animate.set_value(0.5), run_time=1.5,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            sigma.animate.set_value(1), run_time=1.25,
            rate_func=rate_functions.smooth
        )
        self.play(Uncreate(curve))
