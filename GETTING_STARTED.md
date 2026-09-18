# Getting Started (for anyone new to the terminal)

You have a folder of Python files. You're not going to open them and copy-paste code anywhere — you just tell Python to run a file, and it does the work and saves a picture. That's the whole workflow.

## 1. Unzip and open the folder

Unzip `chern-calculator.zip` if you haven't already. Then open VS Code, go to **File > Open Folder**, and select the `chern-calculator` folder (not a file inside it, the whole folder).

## 2. Open a terminal

**Terminal > New Terminal** in the top menu. A panel opens at the bottom with a blinking cursor — that's where you type commands and press Enter.

## 3. Check Python works

```
python --version
```

If that gives you a version number, you're set — use `python` for everything below. If it says "command not found," try `python3 --version` instead, and use `python3` everywhere below if that works. If neither works, install Python from python.org (on Windows, check the box that adds it to PATH during install), then reopen VS Code.

## 4. Install the two packages you need

```
python -m pip install numpy matplotlib
```

One-time setup. Using `python -m pip` instead of just `pip` avoids a common issue where the install goes to a different Python than the one running your scripts.

## 5. Run the scripts

```
python test_models.py            # should print all PASS
python run_bands.py              # saves figures/bands.png
python run_chern.py              # prints Chern numbers, saves the curvature map
python run_phase_diagram.py      # saves the phase diagram
python explore.py                # interactive sliders — opens a window
```

Same idea for the Haldane model files (`run_haldane.py`, `explore_haldane.py`, `run_edge_states.py`).

Open anything in `figures/` from the VS Code sidebar to view it.

## If something breaks

- **`command not found`** — you probably need `python3` instead of `python` (see step 3).
- **`No such file or directory`** — your terminal isn't in the project folder. Run `ls` to check what's there; if `run_bands.py` isn't listed, redo step 1.
- **`ModuleNotFoundError: No module named 'numpy'`** — the install in step 4 went to a different Python. Re-run it with whichever word (`python` or `python3`) you're using to run scripts.
- **No window pops up for `explore.py`** — that's a display/environment issue, not a bug. The other scripts still save their output as images either way.
- Anything else: paste the full terminal output (the command plus everything below it) somewhere you can get help with it — the exact text is what actually pins down the problem.
