"""
Proyecto TESS
Implementación basada en los tutoriales publicados en u-cursos.
Escrito en Python 3, para obtener graficos y resultados, desde una consola GNU/LINUX ejecutar: 
	$ python3 tess_project.py

"""

import numpy as npy
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle
from scipy import stats

def data_handling():
	"""
	Funcion principal del programa, estudia curva de luz del sistema HD74423 contenida en el archivo .fits descargado desde el catalogo MAST.

	"""

	tess_archive = 'tess2019112060037-s0011-0000000355151781-0143-s_lc.fits'
	hdu_list = fits.open(tess_archive)
	light_curve = hdu_list[1].data
	hdu_list.close()

	time = light_curve['TIME']
	flux = light_curve['SAP_FLUX']
	flux_error = light_curve['SAP_FLUX_ERR']
	plot_components = npy.array([time, flux, flux_error])

	normalized_data = data_clean(plot_components)
	npy.sort(normalized_data[1])
	curve_graph_not_clean(plot_components)
	curve_graph(normalized_data)

	frequency, power = LombScargle(normalized_data[0], normalized_data[1], normalized_data[2]).autopower()
	ls_graph(frequency, power)
	npy.argsort(power)[::-1][:100]
	peaks = npy.array([166, 94768, 189703]) 
	best_frequencies = frequency[peaks]
	best_periods = 1 / best_frequencies
	print('Periodos encontrados:', best_periods)
	data_phasing(normalized_data[0], normalized_data[1], best_periods[0], best_periods[1], best_periods[2])

def data_clean(plot_components):
	"""
	Funcion encargada de procesar los componentes a graficar, devolviendo los datos normalizados y sin valores fuera de rango.
	
	"""

	npy.savetxt('light_curve.dat', plot_components)
	clean_components = npy.array([element[~npy.isnan(plot_components[1])] for element in plot_components])
	npy.savetxt('light_curve_filtered.dat', clean_components)
	aux_array = whitout_outliers(clean_components)
	no_outlier_components = npy.array([element[~npy.isnan(aux_array[1])] for element in aux_array])
	npy.savetxt('light_curve_no_outlier.dat', no_outlier_components)
	max_element = npy.amax(no_outlier_components[1])
	no_outlier_components[1] = no_outlier_components[1] / max_element
	npy.savetxt('normalized_data.dat', no_outlier_components)
	return no_outlier_components

def data_phasing(x_1, x_2, phase_1, phase_2, phase_3, *yerr):    
	"""
	Funcion encargada de realizar el ajuste de fase para cierta data y luego graficar el resultado de esta operación.
	
	"""

	phased_t_1 = x_1 % phase_1
	order_1 = npy.argsort(phased_t_1)
	phased_t_2 = x_1 % phase_2
	order_2 = npy.argsort(phased_t_2)
	phased_t_3 = x_1 % phase_3
	order_3 = npy.argsort(phased_t_3)


	if yerr:
		fig, axs = plt.subplots(3, sharex=False, sharey=False)
		fig.suptitle('Ajuste de fase para cada uno de los peaks')
		axs[0].plot(phased_t_1[order_1], x_2[order_1], yerr[0][order_1], 'bo')
		axs[0].legend(['Primer peak'],loc = 1)
		axs[1].plot(phased_t_2[order_2], x_2[order_2], yerr[0][order_2], 'bo')
		axs[1].legend(['Segundo peak'],loc = 1)
		axs[2].plot(phased_t_3[order_3], x_2[order_3], yerr[0][order_3], 'bo')
		axs[2].legend(['Tercer peak'],loc = 1)
		plt.show()
		fig.savefig("exit_3.eps")
	else:
		fig, axs = plt.subplots(3, sharex=False, sharey=False)
		fig.suptitle('Ajuste de fase para cada uno de los peaks')
		axs[0].set_ylabel('magnitud', fontsize = 9)
		axs[0].set_xlabel('fase', fontsize = 9)
		axs[0].plot(phased_t_1[order_1], x_2[order_1], 'bo')
		axs[0].legend(['Primer peak'],loc = 2)
		axs[1].set_ylabel('magnitud', fontsize = 9)
		axs[1].set_xlabel('fase', fontsize = 9)
		axs[1].plot(phased_t_2[order_2], x_2[order_2], 'bo')
		axs[1].legend(['Segundo peak'],loc = 2)
		axs[2].set_ylabel('magnitud', fontsize = 9)
		axs[2].set_xlabel('fase', fontsize = 9)
		axs[2].plot(phased_t_3[order_3], x_2[order_3], 'bo')
		axs[2].legend(['Tercer peak'],loc = 2)
		plt.show()
		fig.savefig("exit_3.eps")

		print('Primer peak: ',phased_t_1[order_1], x_2[order_1])
		print('Segundo peak: ',phased_t_2[order_2], x_2[order_2])
		print('Tercer peak: ',phased_t_3[order_3], x_2[order_3])


def curve_graph(curve, *opt_arg):
	"""
	Función encargada de plotear la curva de luz del sistema, recibe un arreglo como primer argumento a graficar, y un segundo argumento opcional que corresponde al arreglo que servira como ajuste para la primera curva que se grafico.

	"""
	
	figure = plt.figure(figsize=(12, 10))
	axes = plt.gca()
	axes.plot(curve[0], curve[1], alpha=0.9, label='curva de luz')
	axes.set_title('Curva de luz para el sistema estelar HD74423', fontsize=24)
	axes.set_xlabel('Tiempo [d]', fontsize = 18)
	axes.set_ylabel('Amplitud [e- / s]', fontsize = 18)

	if opt_arg:
		popt = opt_arg[0]
		dib_fit = axes.plot(curve[0], sin_fit(curve[0], popt[0], popt[1], popt[2]), alpha=0.9, label='ajuste')

	axes.legend(loc = 4)
	plt.show()
	figure.savefig("exit_1.eps")

def curve_graph_not_clean(curve, *opt_arg):
	"""
	Función encargada de plotear la curva de luz del sistema, recibe un arreglo como primer argumento a graficar, y un segundo argumento opcional que corresponde al arreglo que servira como ajuste para la primera curva que se grafico.

	"""
	
	figure = plt.figure(figsize=(12, 10))
	axes = plt.gca()
	axes.plot(curve[0], curve[1], alpha=0.9, label='curva de luz')
	axes.set_title('Curva de luz para el sistema estelar HD74423 con valores fuera de rango', fontsize=24)
	axes.set_xlabel('Tiempo [d]', fontsize = 18)
	axes.set_ylabel('Amplitud [e- / s]', fontsize = 18)

	if opt_arg:
		popt = opt_arg[0]
		dib_fit = axes.plot(curve[0], sin_fit(curve[0], popt[0], popt[1], popt[2]), alpha=0.9, label='ajuste')

	axes.legend(loc = 4)
	plt.show()
	figure.savefig("exit_0.eps")

def ls_graph(frequency, power):
	"""
	Función encargada de graficar el periodograma de lomb scargle encontrado, se hace 'zoom' en los 3 peaks encontrados.

	"""

	fig, axs = plt.subplots(2, 2)
	axs[0, 0].set_xlabel('Frecuencia [$d^{-1}$]', fontsize = 9)
	axs[0, 0].set_ylabel('Potencia de Lomb-Scargle', fontsize = 9)
	axs[0, 0].plot(frequency, power)
	axs[0, 0].set_title("Periodograma del espectro")
	axs[1, 0].set_xlabel('Frecuencia [$d^{-1}$]', fontsize = 9)
	axs[1, 0].set_ylabel('Potencia de Lomb-Scargle', fontsize = 9)
	axs[1, 0].plot(frequency[0:500], power[0:500])
	axs[1, 0].set_title("Enfoque al primer peak")
	axs[0, 1].set_xlabel('Frecuencia [$d^{-1}$]', fontsize = 9)
	axs[0, 1].set_ylabel('Potencia de Lomb-Scargle', fontsize = 9)
	axs[0, 1].plot(frequency[94700:95200], power[94700:95200])
	axs[0, 1].set_title("Enfoque al segundo peak")
	axs[1, 1].set_xlabel('Frecuencia [$d^{-1}$]', fontsize = 9)
	axs[1, 1].set_ylabel('Potencia de Lomb-Scargle', fontsize = 9)
	axs[1, 1].plot(frequency[189600:190100], power[189600:190100]) 
	axs[1, 1].set_title("Enfoque al tercer peak")
	fig.tight_layout()
	plt.show()
	fig.savefig("exit_2.eps")

	'''
	fig2, axs2 = plt.subplots(3, sharex = True)
	fig2.suptitle('Zoom en peaks de frecuencias para la curva de luz')
	axs2[0].plot(frequency[161:171], power[161:171])
	axs2[1].plot(frequency[94762:94772], power[94762:94772]) 
	axs2[2].plot(frequency[189698:189708], power[189698:189708]) 
	plt.show()
	'''

def whitout_outliers(points, limit = 61500):
	"""
	Elimina aquellos valores que se encuentran fuera de rango para poder hacer un analisis mas focalizado.

	"""

	for i in range(points.shape[1]):

		if points[1][i] < limit:
			points[1][i] = npy.nan

	no_outliers = npy.array([element[~npy.isnan(points[1])] for element in points])

	return no_outliers

def sin_fit(x, w, y, constant):
	"""
	Funcion que por medio de la funcion seno se encarga de ajustar la curva a graficar.

	"""

	return w * npy.sin((2 * npy.pi) * x / y) + constant

if __name__ == '__main__':
	"""
	Se inicia al ejecutar el programa, hace un llamado a la funcion data_handling.
	
	"""

	data_handling()
