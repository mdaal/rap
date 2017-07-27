from setuptools import setup

def readme():
	with open('README.md') as f:
		return f.read()

setup(name='rap',
	version='0.1',
	description='KAM',
	url='http://github.com/storborg/funniest',
	author='mdaal',
	author_email='',
	license='MIT',
	packages=['rap'],
	install_requires=['numpy','scipy' ],
	zip_safe=False)