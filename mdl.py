#------------------------------------------------------#
# magnetic dead layer (MDL) fit plot v0.1		  -----#
# author: tatsunootoshigo, 7475un00705hi90@gmail.com   #
#------------------------------------------------------#

#imports
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes, mark_inset

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams["font.serif"] = "STIX"
mpl.rcParams["mathtext.fontset"] = "stix"

# script version
version = '0.1'
version_name = 'mdl_fit_plot_' + version + '.py'

# parameter adjustment step
valstep_R0 = 1e-2

dataset_name = 'WnFe-SiAlOx'
plot_title = 'MDL_' + dataset_name

sample_R0 = np.array([877.0, 632.0, 388.0, 174.0, 139.0, 83.0])
sample_rho = np.array([88.0, 105.0, 129.0, 174.0, 231.0, 194.0])
FM_thickness = np.array([3, 5, 10, 30, 50, 70])
MR_ratio = ([0.82, 0.609, 0.371, 0.183, 0.109, 0.118],[0.026, 0.035, 0.053, 0.016, 0.138, 0.038],[0.793, 0.584, 0.377, 0.176, 0.238, 0.142])
Delta_R_R0 = ([0.456, 0.235, 0.092, 0.020, 0.01, 0.006],[0.015, 0.014, 0.013, 0.002, 0.012, 0.002],[0.442, 0.235, 0.084, 0.02, 0.021, 0.008]) 
FM_Ms = ([0.7109771, 0.7619016, 0.7858823, 0.703598, 0.718513, 0.6418686])
FM_thickness_range = np.arange(-10.0, 75.0, 1.0)

FM_thickness_cut = np.array([3, 5, 10])
FM_Ms_cut = ([0.7109771, 0.7619016, 0.7858823])
FM_thickness_range_10 = np.arange(-10.0, 15.0, 1.0)

xmin = 0
xmax = 80

ymin = 0
ymax = 50

xprec = 0
yprec = 1

desc_x = 0.1
desc_y = 0.1

out_pdf = 'MDL_fit_plot_' + dataset_name + '.pdf'
out1_pdf = 'MDL_fit_plot_pub_' + dataset_name + '.pdf'
out_svg = 'MDL_fit_plot_' + dataset_name + '.svg'
out1_svg = 'MDL_fit_plot_pub_' + dataset_name + '.svg'

# plot label defs
axis_label_theta = r'$\theta\, / \, \circ$'
axis_label_nm = r'$d_F\;/\;nm$'
axis_label_ohm = r'$R\;/\;\Omega$'
axis_label_msd = r'$\mu_0M_s\cdot d_F\;/\;T\cdot nm$'
axis_label_delta_Ryz = r'$\Delta R_{SMR}\;/\;\%$'

label_xy = r'$xy$'
label_xy_fit = r'$fit$'
label_mdl = r'$MDL$'

# formatting the plot axes
def custom_axis_formater(custom_title, custom_x_label, custom_y_label, xmin, xmax, ymin, ymax, xprec, yprec):
	
	# get axes and tick from plot 
	ax = plt.gca()

	# set the number of major and minor ticks for x,y axes
	# prune='lower' --> remove lowest tick label from x axis
	xmajorLocator = MaxNLocator(12, prune='lower') 
	xmajorFormatter = FormatStrFormatter('%.'+ np.str(xprec) + 'f')
	xminorLocator = MaxNLocator(24) 
	
	ymajorLocator = MaxNLocator(12) 
	ymajorFormatter = FormatStrFormatter('%.'+ np.str(yprec) + 'f')
	yminorLocator = MaxNLocator(24)

	ax.xaxis.set_major_locator(xmajorLocator)
	ax.yaxis.set_major_locator(ymajorLocator)

	ax.xaxis.set_major_formatter(xmajorFormatter)
	ax.yaxis.set_major_formatter(ymajorFormatter)

	# for the minor ticks, use no labels; default NullFormatter
	ax.xaxis.set_minor_locator(xminorLocator)
	ax.yaxis.set_minor_locator(yminorLocator)
	
	# format major and minor ticks width, length, direction 
	ax.tick_params(which='both', width=1, direction='in', labelsize=20)
	ax.tick_params(which='major', length=6)
	ax.tick_params(which='minor', length=4)

	# set axes thickness
	ax.spines['top'].set_linewidth(1.5)
	ax.spines['bottom'].set_linewidth(1.5)
	ax.spines['right'].set_linewidth(1.5)
	ax.spines['left'].set_linewidth(1.5)

	# grid and axes are drawn below the data plot
	ax.set_axisbelow(True)

	# add x,y grids to plot area
	ax.xaxis.grid(True, zorder=0, color='whitesmoke', linestyle='-', linewidth=1)
	ax.yaxis.grid(True, zorder=0, color='whitesmoke', linestyle='-', linewidth=1)

	# set axis labels
	ax.set_xlabel(custom_x_label, fontsize=20)
	ax.set_ylabel(custom_y_label, fontsize=20)

	# set plot title
	#ax.set_title(custom_title, loc='right', fontsize=14)

	return;
def custom_axis_formater_inset(custom_title, custom_x_label, custom_y_label, xmin, xmax, ymin, ymax, xprec, yprec):
	
	# get axes and tick from plot 
	ax = plt.gca()
	# set the number of major and minor bins for x,y axes
	# prune='lower' --> remove lowest tick label from x axis
	xmajorLocator = MaxNLocator(6, prune='lower') 
	xmajorFormatter = FormatStrFormatter('%.'+ np.str(xprec) + 'f')
	xminorLocator = MaxNLocator(12) 
	
	ymajorLocator = MaxNLocator(6) 
	ymajorFormatter = FormatStrFormatter('%.'+ np.str(yprec) + 'f')
	yminorLocator = MaxNLocator(12)
	
	# format major and minor ticks width, length, direction 
	ax.tick_params(which='both', width=1, direction='in', labelsize=20)
	ax.tick_params(which='major', length=6)
	ax.tick_params(which='minor', length=4)

	# set axes thickness
	ax.spines['top'].set_linewidth(1.5)
	ax.spines['bottom'].set_linewidth(1.5)
	ax.spines['right'].set_linewidth(1.5)
	ax.spines['left'].set_linewidth(1.5)

	ax.xaxis.set_major_locator(xmajorLocator)
	ax.yaxis.set_major_locator(ymajorLocator)

	ax.xaxis.set_major_formatter(xmajorFormatter)
	ax.yaxis.set_major_formatter(ymajorFormatter)

	# for the minor ticks, use no labels; default NullFormatter
	ax.xaxis.set_minor_locator(xminorLocator)
	ax.yaxis.set_minor_locator(yminorLocator)

	# grid and axes are drawn below the data plot
	ax.set_axisbelow(True)

	# convert x axis units to radians
	#ax.convert_xunits(radians)

	# add x,y grids to plot area
	ax.xaxis.grid(True, zorder=0, color='gainsboro', linestyle='-', linewidth=1)
	ax.yaxis.grid(True, zorder=0, color='gainsboro', linestyle='-', linewidth=1)

	# set axis labels
	#ax.set_xlabel(custom_x_label, fontsize=20)
	#ax.set_ylabel(custom_y_label, fontsize=20)

	# set plot title
	#ax.set_title(custom_title, loc='right', fontsize=12)

	return;

slope_mdl, intercept_mdl, r_value_mdl, p_value_mdl, std_err_mdl = stats.linregress(FM_thickness[0:4], FM_Ms[0:4]*FM_thickness[0:4])
#print('=================================')
#	print('Ms12: ', intercept_q12, 'Ms34: ', intercept_q34, 'Ms_mean: ', 0.5*(np.absolute(intercept_q12) + np.absolute(intercept_q34)))

print('---------------------------------')
print('slope_mdl: ', slope_mdl)
print('intercept_mdl: ', intercept_mdl)
print('r_value_mdl: ', r_value_mdl)
print('std_err_mdl:', std_err_mdl)
print('---------------------------------')

slope_mdl_cut, intercept_mdl_cut, r_value_mdl_cut, p_value_mdl_cut, std_err_mdl_cut = stats.linregress(FM_thickness[0:2], FM_Ms[0:2]*FM_thickness[0:2])
#print('=================================')
#	print('Ms12: ', intercept_q12, 'Ms34: ', intercept_q34, 'Ms_mean: ', 0.5*(np.absolute(intercept_q12) + np.absolute(intercept_q34)))

print('---------------------------------')
print('slope_mdl: ', slope_mdl_cut)
print('intercept_mdl: ', intercept_mdl_cut)
print('r_value_mdl: ', r_value_mdl_cut)
print('std_err_mdl:', std_err_mdl_cut)
print('---------------------------------')

# create the fig and plot
fig, ax = plt.subplots(figsize=(9, 9), dpi=72)
spec = gridspec.GridSpec(ncols=1, nrows=1)
fig.canvas.set_window_title('mdl_fit_plot_' + version_name) 
plt.figtext(0.80, 0.98, version_name, size=12)

#tx2, = plt.plot(vsm_data[0], fit_params[0]*vsm_data[0] + fit_params[2], 'b--')

tx1, = plt.plot(FM_thickness, FM_Ms*FM_thickness,'ko', mfc='lightgray', markersize=6, label=label_mdl)
tx2, = plt.plot(FM_thickness_range, FM_thickness_range*slope_mdl + intercept_mdl,'k-', mfc='lightgray', markersize=6, label=label_mdl)
# display the legend for the defined labels
plt.legend([tx1, tx2], [label_mdl, 'fit'], loc='upper left', fontsize=14 , frameon=True)
plt.figtext(0.05, 0.92, r'$MDL:$ ' + np.str(np.round(-1.0*intercept_mdl / slope_mdl,3)) + r'$\;nm$' + '	a: ' + np.str(slope_mdl) + '	b: ' + np.str(intercept_mdl) + '\n' + r'$R^2:$' + np.str(r_value_mdl*r_value_mdl) + '	S: ' + np.str(std_err_mdl), size=14)

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

custom_axis_formater(plot_title, axis_label_nm, axis_label_msd, xmin, xmax, ymin, ymax, xprec, yprec)
#-----------------------------------------------------------------------------------------------------

fig1 = plt.figure(figsize=(9, 9), dpi=72)
fig1.canvas.set_window_title(version_name)
spec1 = gridspec.GridSpec(ncols=1, nrows=1)
xy1 = fig1.add_subplot(spec1[0,0])
tx3, = plt.plot(FM_thickness_range, FM_thickness_range*slope_mdl_cut + intercept_mdl_cut,'r-', mfc='lavenderblush', markersize=6, label=label_mdl)
tx4, = plt.plot(FM_thickness, FM_Ms*FM_thickness,'ro', mfc='lavenderblush', markersize=6, label=label_mdl)


# display the legend for the defined labels
plt.legend([tx4, tx3], ['expt', 'fit'], loc='upper left', fontsize=20 , frameon=True)

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

custom_axis_formater(plot_title, axis_label_nm, axis_label_msd, xmin, xmax, ymin, ymax, xprec, yprec)
#xy1.xaxis.grid(True, zorder=0, color='gainsboro', linestyle='-', linewidth=0)
#xy1.yaxis.grid(True, zorder=0, color='gainsboro', linestyle='-', linewidth=0)

fig1.tight_layout(pad=1.0, w_pad=6.0, h_pad=1.0)
plt.subplots_adjust(left=0.1, bottom=0.07, wspace=0.0, hspace=0.0)

axins = inset_axes(xy1,
   						width="40%", # width = 30% of parent_bbox
                        height="40%", # height : 1 inch
                        loc="lower right",
                        borderpad=5)
axins.plot(FM_thickness_range, FM_thickness_range*slope_mdl_cut + intercept_mdl_cut, 'r-', lw=1.5, markersize=8, label=label_xy_fit)
axins.plot(FM_thickness, FM_Ms*FM_thickness, 'ro', mfc='lavenderblush', markersize=8, label=label_xy)
custom_axis_formater_inset(plot_title, axis_label_nm, axis_label_delta_Ryz, xmin, xmax, ymin, ymax, xprec, yprec)

plt.xlim(0, 12)
plt.ylim(0, 12)

#-----------------------------------------------------------------------------------------------------

# create the fig and plot
fig2, ax2 = plt.subplots(figsize=(9, 9), dpi=72)
spec2 = gridspec.GridSpec(ncols=1, nrows=1)
fig2.canvas.set_window_title('mdl_fit_plot_' + version_name) 
plt.figtext(0.80, 0.98, version_name, size=12)

#tx2, = plt.plot(vsm_data[0], fit_params[0]*vsm_data[0] + fit_params[2], 'b--')

tx5, = plt.plot(FM_thickness[0:3], FM_Ms[0:3]*FM_thickness[0:3],'ko', mfc='lightgray', markersize=6, label=label_mdl)
tx6, = plt.plot(FM_thickness_range_10, FM_thickness_range_10*slope_mdl_cut + intercept_mdl_cut,'k-', mfc='lightgray', markersize=6, label=label_mdl)
# display the legend for the defined labels
plt.legend([tx5, tx6], [label_mdl, 'fit'], loc='upper left', fontsize=14 , frameon=True)
plt.figtext(0.05, 0.92, r'$MDL:$ ' + np.str(np.round(-1.0*intercept_mdl_cut / slope_mdl_cut,3)) + r'$\;nm$' + '	a: ' + np.str(slope_mdl_cut) + '	b: ' + np.str(intercept_mdl_cut) + '\n' + r'$R^2:$' + np.str(r_value_mdl_cut*r_value_mdl_cut) + '	S: ' + np.str(std_err_mdl_cut), size=14)

plt.xlim(0, 12)
plt.ylim(0, 12)

custom_axis_formater(plot_title, axis_label_nm, axis_label_msd, xmin, xmax, ymin, ymax, xprec, yprec)
#-----------------------------------------------------------------------------------------------------

pp = PdfPages(out_pdf)
pp.savefig(fig)
pp.close()
pp1 = PdfPages(out1_pdf)
pp1.savefig(fig1)
pp1.close()

# save as .svg too
fig = plt.savefig(out_svg)
fig1 = plt.savefig(out1_svg)

plt.show()