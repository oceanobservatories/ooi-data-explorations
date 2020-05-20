from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name = 'ooi_data_explorations',
    version = '0.0.1',
    description = (
        'Collection of python processing modules for requesting data '
        'from the OOI M2M system'
    ),
    long_description = readme(),
    classifiers = [
        'Development Status:: 1 - Planning',
        'License :: OSI Approved:: MIT License',
        'Programming Language:: Python:: 3:: ONLY',
        'Operating System:: OS Independent',
        'Intended Audience:: Science/Research',
        'Topic:: Scientific/Engineering'
    ],
    keywords = [
        'Ocean Observatories Initiative', 'Regional Cabled Array',
        'Coastal Endurance Array', 'Coastal Pioneer Array',
        'Global Papa Array', 'Global Irminger Array'
    ],
    url = 'https://github.com/oceanobservatories/ooi-data-explorations/python/',
    author = 'Christopher Wingard',
    author_email = 'chris.wingard@oregonstate.edu',
    license = 'MIT',
    packages = find_packages(),
    install_requires = [
        'xarray',
        'munch',
        'tqdm',
        'urllib3',
        'numpy',
        'pandas',
        'gsw',
        'requests',
        'beautifulsoup4',
        'PyYAML'
    ],
    include_package_data=True,
    zip_safe=False
)
