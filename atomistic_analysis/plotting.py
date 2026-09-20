"""LovelyPlots presentation shared by scientific exports and the Anvil browser."""

from contextlib import contextmanager
from importlib.metadata import version, distribution


def style_path():
    # LovelyPlots 1.0.2 registration calls a private API removed in Matplotlib 3.11.
    # Read its installed stylesheet using Matplotlib's public path API instead.
    return str(
        distribution("lovelyplots").locate_file("lovelyplots/styles/ipynb.mplstyle")
    )


def browser_style():
    from matplotlib import rc_params_from_file

    values = rc_params_from_file(style_path(), use_default_template=False)
    return dict(
        provider="LovelyPlots",
        version=version("lovelyplots"),
        style="ipynb",
        colors=values["axes.prop_cycle"].by_key()["color"],
        font_size=values["font.size"],
        font_family="DejaVu Sans, sans-serif",
        line_width=2.5,
        grid=values["axes.grid"],
    )


@contextmanager
def plot_style():
    import matplotlib.pyplot as plt

    with (
        plt.style.context(style_path()),
        plt.rc_context({"lines.linewidth": 2.5, "axes.labelweight": "bold"}),
    ):
        yield
