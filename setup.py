from setuptools import setup, find_packages

setup(
    name='run_dbcan_new',
    version='5.0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'psutil',
        'tqdm',
        'pandas',
        'requests',
        'matplotlib',
        'openpyxl',
        'pyhmmer',
        'biopython==1.84',
        'bcbio-gff',
        'pyrodigal'        
    ],
    entry_points={
        'console_scripts': [
            'run_dbcan=dbcan.run_dbcan:main'
        ]
    },
    author='Xinpeng Zhang',
    author_email='xzhang55@huskers.unl.edu',
    description='Update version for dbCAN',
#    long_description=open('README.md').read(),
#    long_description_content_type='text/markdown',
    url='https://github.com/Xinpeng021001/dbCAN-xinpeng',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)
