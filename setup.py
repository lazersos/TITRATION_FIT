from distutils.core import setup

setup(name='TITRATION_FIT',
	version = '1.0',
	description = 'Library for fitting titration data.',
	long_description = 'This software package contains python '+ \
		'software fitting titration data. ',
	author = 'Samuel A. Lazerson',
	author_email = 'lazersos@gmail.com',
	url = 'http://github.com:lazersos/TITRATION_FIT',
	py_modules=['titration_fit','scipy','matplotlib'],
	scripts = ['titration_fit.py'],
	install_requires=['numpy','scipy','matplotlib']
	)
