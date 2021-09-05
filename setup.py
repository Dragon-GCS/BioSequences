from setuptools import setup, Extension

module1 = Extension(
    'bioseq.algorithm',
    sources=['bioseq/algorithm/algorithmmodule.c', 'bioseq/algorithm/algorithm.c'],
)
# python setup.py bdist_ext --inplace 测试时编译python扩展模块后可正常导入使用

with open("README.md", "r", encoding = "utf-8") as f:
    long_description = f.read()

setup(
        name = 'BioSequences',
        version = '1.0.3',
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
        packages = ['bioseq'],     # or find_packages(exclude=["*.tests", "*.tests.*"...])
        ext_modules=[module1],
)

# reference https://docs.python.org/zh-cn/3.8/distutils/apiref.html#distutils.core.setup
# reference https://packaging.python.org/tutorials/packaging-projects/#packaging-your-project
# 发布方法
# pip install wheel, twine, windows需要vc++14.0以上， linux需要安装python-dev
# python setup.py sdist bdist_wheel 创建源码包和二进制包 linux上添加[--plat-name manylinux2014-x86_64 --universal]
# python -m twine upload dist/*
# 输入用户名和密码