"""plotHeatmap with a layout fix, ported onto deeptools 3.5.6.

lib/deeptools_patched/ vendors a whole deeptools 3.3.2 tree, and that tree no
longer imports: deeptools/cm.py calls matplotlib.cm.register_cmap, which
matplotlib dropped in 3.9. The changes it carries are small, so they are
monkey-patched onto the installed 3.5.6 here instead of maintaining a second
copy of deeptools.

Three patches:

1. prepare_layout - the layout change. The middle region group (Shared)
   gets 1.3x its proportional height, and the summary-plot spacer is
   heatmapwidth/4 instead of /8.

2. --hideProfileGroups LABEL[,LABEL...] - keeps a region group's heatmap rows
   but leaves its line out of the summary profile. With three groups the Shared
   curve mostly traces the envelope of the other two and flattens the contrast
   between them. Implemented by wrapping plot_single, which addProfilePlot
   calls once per group; matplotlib builds the legend from the labelled artists
   on the axis, so a group that is never plotted also leaves the legend.

3. --profileHeight CM - sets the summary-plot height directly. deeptools ties
   it to heatmapwidth, so widening the panels to fit x labels also inflates the
   profile and squashes the heatmaps. Defaults to heatmapwidth, the old
   behaviour, when the flag is absent.

4. --profileSpacer CM - the gap between the summary plot and the heatmaps.
   The patch widened it from heatmapwidth/8 to /4; tying it to the width means it
   grows whenever the panels do. Defaults to heatmapwidth/4 when absent.

5. --heatmapTickSize PT - the heatmap x tick label size. deeptools hardcodes
   size=8 there (plotHeatmap.py line 714), so the local matplotlibrc governs
   the profile ticks but not the heatmap ticks, and narrowing the panels makes
   the heatmap labels collide while the profile ones still fit. Only calls that
   pass an explicit size are touched, which is that one call.

6. --panelWspace FRACTION / --panelWpad INCHES - horizontal room between the
   sample columns. deeptools builds the figure with constrained_layout, so the
   gap comes from that engine (wspace 0.02, w_pad 0.04167 by default) and not
   from the GridSpec, which is why width_ratios cannot close it.

7. --legendOutside - moves the profile legend out of the plot area to the right
   of the last profile panel, vertically centred. Every deeptools
   --legendLocation choice is inside the axes, where the legend sits on top of
   the curves. deeptools makes exactly one legend call, on ax_list[-1], so this
   wraps Axes.legend and disarms itself after that call. savefig already uses
   bbox_inches='tight', so the legend is not clipped.
"""
import sys

import numpy as np
from matplotlib import gridspec
from matplotlib.axes import Axes

from deeptools import plotHeatmap


def _take_flags(argv):
    """Pull our own flags out of argv before deeptools parses it."""
    hide = set()
    outside = False
    prof = None
    spacer = None
    ticksize = None
    wspace = None
    wpad = None
    out = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--hideProfileGroups':
            if i + 1 >= len(argv):
                sys.exit('--hideProfileGroups needs a value')
            hide.update(x for x in argv[i + 1].split(',') if x)
            i += 2
            continue
        if a.startswith('--hideProfileGroups='):
            hide.update(x for x in a.split('=', 1)[1].split(',') if x)
            i += 1
            continue
        if a == '--profileHeight':
            if i + 1 >= len(argv):
                sys.exit('--profileHeight needs a value')
            prof = float(argv[i + 1])
            i += 2
            continue
        if a.startswith('--profileHeight='):
            prof = float(a.split('=', 1)[1])
            i += 1
            continue
        if a == '--profileSpacer':
            if i + 1 >= len(argv):
                sys.exit('--profileSpacer needs a value')
            spacer = float(argv[i + 1])
            i += 2
            continue
        if a.startswith('--profileSpacer='):
            spacer = float(a.split('=', 1)[1])
            i += 1
            continue
        if a == '--heatmapTickSize':
            if i + 1 >= len(argv):
                sys.exit('--heatmapTickSize needs a value')
            ticksize = float(argv[i + 1])
            i += 2
            continue
        if a.startswith('--heatmapTickSize='):
            ticksize = float(a.split('=', 1)[1])
            i += 1
            continue
        if a == '--panelWspace':
            if i + 1 >= len(argv):
                sys.exit('--panelWspace needs a value')
            wspace = float(argv[i + 1])
            i += 2
            continue
        if a.startswith('--panelWspace='):
            wspace = float(a.split('=', 1)[1])
            i += 1
            continue
        if a == '--panelWpad':
            if i + 1 >= len(argv):
                sys.exit('--panelWpad needs a value')
            wpad = float(argv[i + 1])
            i += 2
            continue
        if a.startswith('--panelWpad='):
            wpad = float(a.split('=', 1)[1])
            i += 1
            continue
        if a == '--legendOutside':
            outside = True
            i += 1
            continue
        out.append(a)
        i += 1
    return hide, outside, prof, spacer, ticksize, wspace, wpad, out


HIDE_GROUPS = set()
LEGEND_OUTSIDE = False
PROFILE_HEIGHT = None
PROFILE_SPACER = None
HEATMAP_TICK_SIZE = None
PANEL_WSPACE = None
PANEL_WPAD = None

_plot_single = plotHeatmap.plot_single
_axes_legend = Axes.legend
_axes_set_xticklabels = Axes.set_xticklabels


def plot_single(ax, ma, average_type, color, label, plot_type='lines'):
    if label in HIDE_GROUPS:
        return ax
    return _plot_single(ax, ma, average_type, color, label, plot_type=plot_type)


def legend(self, *args, **kwargs):
    global LEGEND_OUTSIDE
    if LEGEND_OUTSIDE:
        LEGEND_OUTSIDE = False          # only the first legend, the profile one
        kwargs['loc'] = 'center left'
        kwargs['bbox_to_anchor'] = (1.02, 0.5)
    return _axes_legend(self, *args, **kwargs)


def set_xticklabels(self, labels, *args, **kwargs):
    if HEATMAP_TICK_SIZE is not None and ('size' in kwargs or 'fontsize' in kwargs):
        kwargs.pop('fontsize', None)
        kwargs['size'] = HEATMAP_TICK_SIZE
    return _axes_set_xticklabels(self, labels, *args, **kwargs)


def prepare_layout(hm_matrix, heatmapsize, showSummaryPlot, showColorbar,
                   perGroup, colorbar_position, fig):
    heatmapwidth, heatmapheight = heatmapsize

    if PANEL_WSPACE is not None or PANEL_WPAD is not None:
        engine = fig.get_layout_engine()
        if engine is not None:
            opts = {}
            if PANEL_WSPACE is not None:
                opts['wspace'] = PANEL_WSPACE
            if PANEL_WPAD is not None:
                opts['w_pad'] = PANEL_WPAD
            try:
                engine.set(**opts)
            except (AttributeError, TypeError):
                pass          # a layout engine that does not take these

    numcols = hm_matrix.get_num_samples()
    numrows = hm_matrix.get_num_groups()
    if perGroup:
        numcols, numrows = numrows, numcols

    if perGroup:
        height_ratio = np.array([np.amax(np.diff(hm_matrix.group_boundaries))] * numrows).astype(float)
    else:
        height_ratio = np.diff(hm_matrix.group_boundaries).astype(float)
        # make the middle region group (Shared) ~1.3x taller
        if len(height_ratio) >= 3:
            height_ratio[len(height_ratio) // 2] *= 1.3
    height_ratio = heatmapheight * (height_ratio / height_ratio.sum())

    height_ratio = height_ratio.tolist()
    width_ratio = [heatmapwidth] * numcols

    if showColorbar:
        if colorbar_position == 'below':
            numrows += 2
            height_ratio += [4 / 2.54]
            height_ratio += [1 / 2.54]
        else:
            numcols += 1
            width_ratio += [1 / 2.54]

    if showSummaryPlot:
        numrows += 2
        sumplot_height = heatmapwidth if PROFILE_HEIGHT is None else PROFILE_HEIGHT
        spacer_height = heatmapwidth / 4 if PROFILE_SPACER is None else PROFILE_SPACER
        height_ratio = np.concatenate([[sumplot_height, spacer_height], height_ratio])

    return gridspec.GridSpec(numrows, numcols, height_ratios=height_ratio,
                             width_ratios=width_ratio, figure=fig)


plotHeatmap.prepare_layout = prepare_layout
plotHeatmap.plot_single = plot_single
Axes.legend = legend
Axes.set_xticklabels = set_xticklabels

if __name__ == "__main__":
    (HIDE_GROUPS, LEGEND_OUTSIDE, PROFILE_HEIGHT, PROFILE_SPACER,
     HEATMAP_TICK_SIZE, PANEL_WSPACE, PANEL_WPAD, sys.argv) = _take_flags(sys.argv)
    plotHeatmap.main()
