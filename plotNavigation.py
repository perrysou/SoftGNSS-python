import matplotlib as mpl
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d

import initSettings

# %% configure matplotlib
mpl.rcdefaults()
# mpl.rcParams['font.sans-serif']
# mpl.rcParams['font.family'] = 'serif'
mpl.rc('savefig', bbox='tight', transparent=False, format='png')
mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
mpl.rc('figure', figsize=[8, 6], autolayout=False, dpi=120)
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', serif='computer modern roman', size=10)
mpl.rc('mathtext', fontset='cm')


# mpl.rc('font', size=16)
# mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

# ./plotNavigation.m


def plotNavigation(navSolutions=None, settings=None, *args, **kwargs):
    # Functions plots variations of coordinates over time and a 3D position
    # plot. It plots receiver coordinates in UTM system or coordinate offsets if
    # the true UTM receiver coordinates are provided.

    # plotNavigation(navSolutions, settings)

    #   Inputs:
    #       navSolutions    - Results from navigation solution function. It
    #                       contains measured pseudoranges and receiver
    #                       coordinates.
    #       settings        - Receiver settings. The true receiver coordinates
    #                       are contained in this structure.

    ## Plot results in the necessary data exists ==============================
    if navSolutions is not None:
        refCoord = initSettings.TruePosition()
        ## If reference position is not provided, then set reference position
        ## to the average postion
        if settings.truePosition.E is None or settings.truePosition.N is None or settings.truePosition.U is None:
            # === Compute mean values ==========================================
            # Remove NaN-s or the output of the function MEAN will be NaN.
            refCoord.E = np.nanmean(navSolutions[0].E)

            refCoord.N = np.nanmean(navSolutions[0].N)

            refCoord.U = np.nanmean(navSolutions[0].U)

            meanLongitude = np.nanmean(navSolutions[0].longitude)

            meanLatitude = np.nanmean(navSolutions[0].latitude)

            refPointLgText = 'Mean Position' + '\\newline Lat: %.5f $^\circ$' % meanLatitude + \
                             '\\newline Lng: %.5f $^\circ$' % meanLongitude + \
                             '\\newline Hgt: %+6.1f' % np.nanmean(navSolutions[0].height)

        else:
            refPointLgText = 'Reference Position'

            refCoord.E = settings.truePosition.E

            refCoord.N = settings.truePosition.N

            refCoord.U = settings.truePosition.U

        figureNumber = 300

        # figure windows, when many figures are closed and reopened. Figures
        # drawn or opened by the user, will not be "overwritten" by this
        # function if the auto numbering is not used.
        # === Select (or create) and clear the figure ==========================
        f = plt.figure(figureNumber)
        f.clf()
        f.set_label('Navigation solutions')
        spec = gs.GridSpec(2, 2)
        h11 = plt.subplot(spec[0:2])

        # the axes3d module is needed for the following line
        dummy = axes3d.Axes3D
        h31 = plt.subplot(spec[2], projection='3d')

        h32 = plt.subplot(spec[3], projection='polar')

        ## Plot all figures =======================================================
        # --- Coordinate differences in UTM system -----------------------------
        h11.plot(navSolutions[0].E - refCoord.E, '-',
                 navSolutions[0].N - refCoord.N, '-',
                 navSolutions[0].U - refCoord.U, '-')
        h11.legend(['E', 'N', 'U'])
        h11.set(title='Coordinates variations in UTM system',
                xlabel='Measurement period: %i ms' % settings.navSolPeriod,
                ylabel='Variations (m)')
        h11.grid()
        h11.axis('tight')
        h31.plot((navSolutions[0].E - refCoord.E).T,
                 (navSolutions[0].N - refCoord.N).T,
                 (navSolutions[0].U - refCoord.U).T, '+')
        h31.hold(True)
        h31.plot([0], [0], [0], 'r+', lw=1.5, ms=10)
        h31.hold(False)
        # h31.viewLim(0,90)
        h31.axis('equal')
        h31.grid(which='minor')
        h31.legend(['Measurements', refPointLgText])
        h31.set(title='Positions in UTM system (3D plot)',
                xlabel='East (m)',
                ylabel='North (m)',
                zlabel='Upping (m)')
        h32.plot(np.deg2rad(navSolutions[0].channel[0].az.T),
                 90 - navSolutions[0].channel[0].el.T)
        [h32.text(x, y, s) for x, y, s in zip(np.deg2rad(navSolutions[0].channel[0].az[:, 0]),
                                              90 - navSolutions[0].channel[0].el[:, 0],
                                              navSolutions[0].channel[0].PRN[:, 0])]
        h32.set_theta_direction(-1)
        h32.set_theta_zero_location('N')
        h32.set_xlim([0, 2 * np.pi])
        h32.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))
        h32.set_rlabel_position(0)
        h32.set_ylim([0, 90])
        h32.set_yticks([0, 15, 30, 45, 60, 75])
        h32.set_yticklabels([90, 75, 60, 45, 30, 15])
        h32.set_title('Sky plot (mean PDOP: %f )' % np.mean(navSolutions[0].DOP[1, :]))
        f.show()
    else:
        print 'plotNavigation: No navigation data to plot.'


if __name__ == '__main__':
    navSolutions = np.load('./navSolutions_python.npy')
    settings = initSettings.Settings()
    plotNavigation(navSolutions, settings)
