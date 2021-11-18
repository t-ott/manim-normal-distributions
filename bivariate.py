from manim import *
import math

# Nice color ramp, cool (purplish) to hot (yellowish) RGB color values
COLOR_RAMP = [
    rgb_to_color([57/255, 0.0, 153/255]),
    rgb_to_color([158/255, 0.0, 89/255]),
    rgb_to_color([1.0, 0.0, 84/255]),
    rgb_to_color([1.0, 84/255, 0.0]),
    rgb_to_color([1.0, 189/255, 0.0])
]


def PDF_bivariate_normal(x_1, x_2, mu_1=0, mu_2=0, sigma_1=1, sigma_2=1, rho=0):
    '''
    General form of probability density function of bivariate normal distribution
    '''
    normalizing_const = 1/(2 * math.pi * sigma_1 * sigma_2 * math.sqrt(1 - rho**2))
    exp_coeff = -(1/(2 * (1 - rho**2)))
    A = ((x_1 - mu_1)/sigma_1)**2
    B = -2 * rho * ((x_1 - mu_1)/sigma_1) * ((x_2 - mu_2)/sigma_2)
    C = ((x_2 - mu_2)/sigma_2)**2

    return normalizing_const * math.exp(exp_coeff*(A + B + C))


class StandardBivariateNormal(ThreeDScene):
    '''
    Plots the surface of the probability density function of the standard
    bivariate normal distribution
    '''

    def construct(self):
        ax = ThreeDAxes(
            x_range = [-4, 4, 1],
            y_range = [-4, 4, 1],
            z_range = [0, 0.2, 0.1]
        )
        x_label = ax.get_x_axis_label(r'x_1')
        y_label = ax.get_y_axis_label(r'x_2', edge=UP, buff=0.2)
        z_label = ax.get_z_axis_label(r'\phi(x_1, x_2)', buff=0.2)
        axis_labels = VGroup(x_label, y_label, z_label)

        # Define Surface using the default values of PDF_bivariate_normal()
        # which represent the standard bivariate normal distribution
        distribution = Surface(
            lambda u, v: ax.c2p(u, v, PDF_bivariate_normal(u, v)),
            resolution=(42, 42),
            u_range=[-3.5, 3.5],
            v_range=[-3.5, 3.5],
            fill_opacity=0.7
        )
        distribution.set_fill_by_value(
            axes = ax,
            # Utilize color ramp colors, higher values are "warmer"
            colors = [(COLOR_RAMP[0], 0),
                      (COLOR_RAMP[1], 0.05),
                      (COLOR_RAMP[2], 0.1),
                      (COLOR_RAMP[3], 0.15),
                      (COLOR_RAMP[4], 0.2)]
        )

        # Set up animation
        self.add(ax, axis_labels)
        self.set_camera_orientation(
            phi=75*DEGREES,
            theta=-70*DEGREES,
            frame_center=[0, 0, 2],
            zoom=0.75)
        # Begin animation
        self.play(Create(distribution))
        self.move_camera(theta=70*DEGREES, run_time=2)
        self.move_camera(theta=-70*DEGREES, run_time=2)
        self.play(Uncreate(distribution))


class AdjustMu(ThreeDScene):
    '''
    Scene plots the surface of the probability density function of the bivariate
    normal distribution, then animates various adjustments to the means (mu) of
    x_1 and x_2
    '''

    def construct(self):
        ax = ThreeDAxes(
            x_range = [-5, 5, 1],
            y_range = [-5, 5, 1],
            z_range = [0, 0.2, 0.1]
        )
        x_label = ax.get_x_axis_label(r'x_1')
        y_label = ax.get_y_axis_label(r'x_2', edge=UP, buff=0.2)
        z_label = ax.get_z_axis_label(r'\phi(x_1, x_2)', buff=0.2)
        axis_labels = VGroup(x_label, y_label, z_label)

        # Initialize ValueTrackers to adjust means
        mu_1 = ValueTracker(0)
        mu_2 = ValueTracker(0)

        mu_1_tex = MathTex(r'\mu_1 =')
        mu_1_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2, include_sign=True)
            .set_value(mu_1.get_value())
            .next_to(mu_1_tex, RIGHT)
        )
        mu_1_text_group = VGroup(mu_1_tex, mu_1_value_text).arrange(RIGHT, buff=0.5)

        mu_2_tex = MathTex(r'\mu_2 =')
        mu_2_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2, include_sign=True)
            .set_value(mu_2.get_value())
            .next_to(mu_2_tex, RIGHT)
        )
        mu_2_text_group = VGroup(mu_2_tex, mu_2_value_text).arrange(RIGHT, buff=0.5)

        text = VGroup(mu_1_text_group, mu_2_text_group).arrange(RIGHT, buff=2).move_to(UP*3.5)
        # Fix position of text so 3D camera movements do not impact position
        self.add_fixed_in_frame_mobjects(text)

        # Define suface function of PDF, always redraw to allow for smooth
        # value adjustments to mu_1 and mu_2
        distribution = always_redraw(
            lambda: Surface(
                lambda u, v: ax.c2p(
                    u, v, PDF_bivariate_normal(
                        u, v, mu_1=mu_1.get_value(), mu_2=mu_2.get_value()
                    )
                ),
                resolution=(42, 42),
                u_range=[-4.5, 4.5],
                v_range=[-4.5, 4.5],
                fill_opacity=0.7
            ).set_fill_by_value(
                axes = ax,
                # Utilize color ramp colors with, higher values are "warmer"
                colors = [(COLOR_RAMP[0], 0),
                          (COLOR_RAMP[1], 0.05),
                          (COLOR_RAMP[2], 0.1),
                          (COLOR_RAMP[3], 0.15),
                          (COLOR_RAMP[4], 0.2)]
            )
        )

        # Set up animation
        self.set_camera_orientation(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 2],
            zoom=0.6
        )
        self.add(ax, axis_labels, text)
        # Begin animation
        self.play(Create(distribution))
        self.play(
            mu_1.animate.set_value(2), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu_1.animate.set_value(0), run_time=1,
            rate_func=rate_functions.smooth
        )
        self.play(
            mu_2.animate.set_value(-2), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu_2.animate.set_value(0), run_time=1,
            rate_func=rate_functions.smooth
        )
        # Top-down view
        self.move_camera(
            theta=-90*DEGREES,
            phi=0,
            frame_center=[0, 0, 0],
            zoom=0.5
        )
        self.play(
            mu_1.animate.set_value(-2),
            mu_2.animate.set_value(2),
            run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu_1.animate.set_value(0),
            mu_2.animate.set_value(0),
            run_time=1,
            rate_fun=rate_functions.smooth
        )
        # Return camera to original position
        self.move_camera(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 2],
            zoom=0.6
        )
        self.play(Uncreate(distribution))
        self.wait()


class AdjustSigma(ThreeDScene):
    '''
    Scene plots the surface of the probability density function of the bivariate
    normal distribution, then animates various adjustments to the standard
    deviations of x_1 and x_2
    '''

    def construct(self):
        resolution_fa = 50
        ax = ThreeDAxes(
            x_range = [-5, 5, 1],
            y_range = [-5, 5, 1],
            z_range = [0, 0.4, 0.1]
        )
        x_label = ax.get_x_axis_label(r'x_1')
        y_label = ax.get_y_axis_label(r'x_2', edge=UP, buff=0.2)
        z_label = ax.get_z_axis_label(r'\phi(x_1, x_2)', buff=0.2)
        axis_labels = VGroup(x_label, y_label, z_label)

        # Initialize ValueTracker objects to adjust standard deviations
        sigma_1 = ValueTracker(1)
        sigma_2 = ValueTracker(1)

        sigma_1_tex = MathTex(r'\sigma_1 =')
        sigma_1_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2, include_sign=True)
            .set_value(sigma_1.get_value())
            .next_to(sigma_1_tex, RIGHT)
        )
        sigma_1_text_group = VGroup(sigma_1_tex, sigma_1_value_text).arrange(RIGHT, buff=0.5)

        sigma_2_tex = MathTex(r'\sigma_2 =')
        sigma_2_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2, include_sign=True)
            .set_value(sigma_2.get_value())
            .next_to(sigma_2_tex, RIGHT)
        )
        sigma_2_text_group = VGroup(sigma_2_tex, sigma_2_value_text).arrange(RIGHT, buff=0.5)

        text = VGroup(sigma_1_text_group, sigma_2_text_group).arrange(RIGHT, buff=2).move_to(UP*3.5)
        # Fix position of text so 3D camera movements do not impact position
        self.add_fixed_in_frame_mobjects(text)

        # Define suface function of PDF, always redraw to allow for smooth
        # value adjustments to sigma_1 and sigma_2
        distribution = always_redraw(
            lambda: Surface(
                lambda u, v: ax.c2p(
                    u, v, PDF_bivariate_normal(
                        u, v, sigma_1=sigma_1.get_value(), sigma_2=sigma_2.get_value()
                    )
                ),
                resolution=(resolution_fa, resolution_fa),
                u_range=[-4.5, 4.5],
                v_range=[-4.5, 4.5],
                fill_opacity=0.7
            ).set_fill_by_value(
                axes = ax,
                # Utilize color ramp colors with, higher values are "warmer"
                colors = [(COLOR_RAMP[0], 0),
                          (COLOR_RAMP[1], 0.05),
                          (COLOR_RAMP[2], 0.1),
                          (COLOR_RAMP[3], 0.15),
                          (COLOR_RAMP[4], 0.2)]
            )
        )

        self.set_camera_orientation(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 3],
            zoom=0.6
        )
        self.add(ax, axis_labels, text)

        self.play(Create(distribution))
        self.play(
            sigma_1.animate.set_value(0.5),
            run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.move_camera(theta=70*DEGREES, run_time=2)
        self.move_camera(theta=-70*DEGREES, run_time=2)
        self.play(
            sigma_1.animate.set_value(1),
            run_time=1,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            sigma_2.animate.set_value(1.5),
            run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.move_camera(theta=70*DEGREES, run_time=2)
        self.move_camera(theta=-70*DEGREES, run_time=2)
        self.play(
            sigma_2.animate.set_value(1),
            run_time=1,
            rate_func=rate_functions.smooth
        )
        self.play(
            sigma_1.animate.set_value(1.5),
            sigma_2.animate.set_value(0.5),
            run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        # Top-down view
        self.move_camera(
            theta=-90*DEGREES,
            phi=0,
            frame_center=[0, 0, 0],
            zoom=0.5
        )
        self.wait()
        self.play(
            sigma_1.animate.set_value(0.5),
            sigma_2.animate.set_value(1.5),
            run_time=1,
            rate_fun=rate_functions.smooth
        )
        self.wait(2)
        self.play(
            sigma_1.animate.set_value(0.75),
            sigma_2.animate.set_value(0.75),
        )
        # Return camera to starting position
        self.move_camera(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 3],
            zoom=0.6
        )
        self.wait()
        self.play(
            sigma_1.animate.set_value(1),
            sigma_2.animate.set_value(1),
        )
        self.play(Uncreate(distribution))
        self.wait()


class AdjustRho(ThreeDScene):
    '''
    Scene plots the surface of the probability density function of the bivariate
    normal distribution, then animates various adjustments to rho, the
    correlation between x_1 and x_2
    '''

    def construct(self):
        resolution_fa = 55
        ax = ThreeDAxes(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            z_range=[0, 0.3, 0.1]
        )
        x_label = ax.get_x_axis_label(r'x_1')
        y_label = ax.get_y_axis_label(r'x_2', edge=UP, buff=0.2)
        z_label = ax.get_z_axis_label(r'\phi(x_1, x_2)', buff=0.2)
        axis_labels = VGroup(x_label, y_label, z_label)

        # Define rho ValueTracker to animate adjustments to rho
        rho = ValueTracker(0)

        rho_tex = MathTex(r'\rho =')
        rho_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2, include_sign=True)
            .set_value(rho.get_value())
            .next_to(rho_tex, RIGHT)
        )
        text = VGroup(rho_tex, rho_value_text).arrange(RIGHT, buff=0.5).move_to(UP*3.5)
        # Fix position of text so 3D camera movements do not impact position
        self.add_fixed_in_frame_mobjects(text)

        # Define suface function of PDF, always redraw to allow for smooth
        # value adjustments to rho
        distribution = always_redraw(
            lambda: Surface(
                lambda u, v: ax.c2p(
                    u, v, PDF_bivariate_normal(u, v, rho=rho.get_value())
                ),
                resolution=(resolution_fa, resolution_fa),
                u_range=[-4.5, 4.5],
                v_range=[-4.5, 4.5],
                fill_opacity=0.7
            ).set_fill_by_value(
                axes = ax,
                # Utilize color ramp colors, higher values are "warmer"
                colors = [(COLOR_RAMP[0], 0),
                          (COLOR_RAMP[1], 0.05),
                          (COLOR_RAMP[2], 0.1),
                          (COLOR_RAMP[3], 0.15),
                          (COLOR_RAMP[4], 0.2)]
            )
        )

        # Initialize camera position and mobjects
        self.set_camera_orientation(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 2.5],
            zoom=0.6
        )
        self.add(text, ax, axis_labels)

        # Begin animation
        self.play(Create(distribution))
        # Animate to positive rho (positive correlation)
        self.play(
            rho.animate.set_value(0.75), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        # Look down on surface
        self.move_camera(
            theta=-90*DEGREES,
            phi=0,
            frame_center=[0, 0, 0],
            zoom=0.5
        )
        self.wait(2)
        # Return to starting camera position
        self.move_camera(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 2.5],
            zoom=0.6
        )
        # Animate to negative rho (negative correlation)
        self.play(
            rho.animate.set_value(-0.75), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        # Look down on surface
        self.move_camera(
            theta=-90*DEGREES,
            phi=0,
            frame_center=[0, 0, 0],
            zoom=0.5
        )
        self.wait(2)
        # Animate to zero rho (no correlation)
        self.play(
            rho.animate.set_value(0), run_time=2,
            rate_func=rate_functions.smooth
        )
        # Return camera to starting position
        self.move_camera(
            theta=-70*DEGREES,
            phi=70*DEGREES,
            frame_center=[0, 0, 2.5],
            zoom=0.6
        )
        self.play(Uncreate(distribution))
        self.wait()
