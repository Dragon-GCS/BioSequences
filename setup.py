from setuptools import setup, Extension

module1 = Extension(
    'biosequence.algorithm',
    sources=['./biosequence/algorithm/algorithmmodule.c', './biosequence/algorithm/algorithm.c'],
)

with open("README.md", "r", encoding = "utf-8") as f:
    long_description = f.read()

setup(
        name = 'BioSequences',
        version = '1.0.1',
        author = 'Dragon',
        author_email = '878173121@qq.com',
        description = 'Tools to analysis biology sequence',
        long_description=long_description,
        long_description_content_type="text/markdown",
        license = 'GPL Licence',
        url = 'https://github.com/Dragon-GCS/BioSequence',
        project_urls={
                "Bug Tracker": "https://github.com/Dragon-GCS/BioSequence/issues",
            },
        install_requires = ["matplotlib >= 3.4.3"],
        classifiers = [
                'Topic :: Scientific/Engineering :: Bio-Informatics',
                "Environment :: Console",
                'License :: OSI Approved :: GNU General Public License (GPL)',
                'Intended Audience :: Science/Research',
                'Operating System :: OS Independent',
                'Natural Language :: Chinese (Simplified)',
                "Natural Language :: English",
                'Programming Language :: Python :: 3.7',
              ],
        keywords=["biology", "analysis"],
        packages = ['biosequence'],     # or find_packages(exclude=["*.tests", "*.tests.*"...])
        ext_modules=[module1],
)

# reference https://docs.python.org/zh-cn/3.8/distutils/apiref.html#distutils.core.setup
# reference https://packaging.python.org/tutorials/packaging-projects/#packaging-your-project
