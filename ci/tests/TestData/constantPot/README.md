# Sample Potential Hook

This is an example of how one hooks a C++ or C potential into python so that it can be directly used by CPotentialLib

The only changes that need to be made are changing 

```lang-c
double constantPot(std::vector<std::vector<double> >, std::vector<std::string>) {
    return 52.0
}
```

in `constantPot.cpp` to the actual potential you have.

After that, changing `constantPot` to the name you'd like for your potential should suffice to allow it to compile

At that point calling `CPotential(<your_module_name>.potential)` should be enough to set up the C-potential hook.