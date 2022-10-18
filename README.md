# Psience [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/McCoyGroup/Binder-McUtils/master?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FMcCoyGroup%252FPsience%26urlpath%3Dlab%252Ftree%252FPsience%252Fbinder%252Findex.ipynb%26branch%3Dmaster)

Psience is a set of core scientific packages written by the McCoy group for the McCoy group to handle interesting scientific problems, like DVR, managing potential and dipole surfaces, VPT2, normal mode analysis, etc.

We're working on [documenting the package](https://mccoygroup.github.io/Psience), but writing good documentation takes more time than writing good code.

### Installation & Requirements

`Psience` is written in pure python and we've worked to try to avoid any major dependencies outside of what comes in `Anaconda` and our `McUtils` package.

The easiest way to install is via `pip`, as

```lang-shell
pip install mccoygroup-psience
```

This should install all dependencies. 
The major requirement is that Python 3.8+ is required due to use of the types module. If installing on Windows, use  the above command in the "Anaconda Prompt" as opposed to the default terminal installed on Windows.
For safety, it is best to install this in a [virtual environment](https://docs.python.org/3.8/tutorial/venv.html), which we can make like

```lang-shell
python3.8 -m pip venv mcenv
```

and activate like

```lang-shell
. mcenv/bin/activate
```

or to use it in a [container](https://www.docker.com/) or [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or some other place where we can control the environment.

It is also possible to install from source like

```lang-shell
git clone https://github.com/McCoyGroup/Psience.git
```

but in this case you will need to make sure the library is on the path yourself and all of the dependencies are installed.

Once the package is installed, you can go ahead and get started in your scripts by importing Psience with the following command.

```lang-shell
import Psience
```

Have fun doing Psience!

### Contributing

If you'd like to help out with this, we'd love contributions.
The easiest way to get started with it is to try it out.
When you find bugs, please [report them](https://github.com/McCoyGroup/Psience/issues/new?title=Bug%20Found:&labels=bug). 
If there are things you'd like added [let us know](https://github.com/McCoyGroup/Psience/issues/new?title=Feature%20Request:&labels=enhancement), and we'll try to help you get the context you need to add them yourself.
One of the biggest places where people can help out, though, is in improving the quality of the documentation.
As you try things out, add them as examples, either to the [main page](https://mccoygroup.github.io/References/Documentation/Psience.html#examples) or to a [child page](https://mccoygroup.github.io/References/Documentation/Psience/Molecools/Molecule/Molecule.html#examples).
You can also edit the docstrings in the code to add context, explanation, argument types, return types, etc.
