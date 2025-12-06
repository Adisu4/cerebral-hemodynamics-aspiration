from setuptools import setup, find_packages

setup(
    name='cerebral-hemodynamics-aspiration',
    version='1.0.0',
    description='Computational model of cerebral hemodynamics with venous aspiration therapy for ICP management',
    author='Adisu Mengesha Assefa',
    author_email='aassefa@unomaha.edu',
    url='https://github.com/adisumengesha/cerebral-hemodynamics-aspiration',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'matplotlib>=3.3.0',
    ],
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
