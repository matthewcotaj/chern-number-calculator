# Getting Started — assuming you've never used a terminal

This is the click-by-click version. It assumes no command-line experience. If you
follow it top to bottom you'll go from the downloaded zip to seeing all three
figures. Total time: about 20 minutes, most of it a one-time setup.

---

## The single most important idea

You have a **folder of Python files**. You do **not** open those files and paste
their contents into the terminal. The files just sit there. You type one short
command that tells Python: "go run that file for me." Python opens the file,
does the work, and saves a picture.

If you ever see the terminal stuck showing just a `>` or the word `quote>` or
`dquote>` on a new line, it means you accidentally pasted code (or an unfinished
command) into it. Fix: click in the terminal and press **Ctrl + C** once. You're
back to normal. Then run a file the correct way (below) instead of pasting into
it.

---

## Step 1 — Put the folder on your computer

1. Find the `chern-calculator.zip` you downloaded (probably in your Downloads
   folder).
2. Unzip it. On Windows: right-click > "Extract All". On Mac: double-click it.
   You now have a normal folder called `chern-calculator`.
3. Remember where it is. Moving it to your Desktop is fine and makes it easy to
   find.

---

## Step 2 — Open the FOLDER in VS Code

This part matters. Open the whole folder, not a single file.

1. Open VS Code.
2. Top menu: **File > Open Folder…**
3. Select the `chern-calculator` folder and confirm.
4. On the left you should now see a list of the files: `README.md`,
   `run_bands.py`, a `qwz` folder, and so on.

Why this matters: when you open the folder, VS Code's terminal will
automatically start *inside* that folder. That saves you from having to learn how
to navigate folders in the terminal.

---

## Step 3 — Open the terminal inside VS Code

1. Top menu: **Terminal > New Terminal**.
2. A panel opens at the bottom. You'll see a line of text ending in a symbol like
   `$` or `>` or `%`, with a blinking cursor. That symbol is the "prompt." It
   means the terminal is ready for a command.
3. Look at the last part of the prompt line — it should mention
   `chern-calculator`. That confirms the terminal is pointing at your project
   folder. Good.

You type a command, press **Enter**, and it runs. That's the whole interaction.

---

## Step 4 — Check that Python is installed

Type this and press Enter:

```
python --version
```

- If it prints something like `Python 3.11.5`, great. Your python word is
  `python`. Use that everywhere below.
- If it says `command not found` or `Python was not found`, try:

```
python3 --version
```

- If *that* prints a version, your python word is `python3`. Use `python3`
  everywhere below instead of `python`.
- If **both** fail, Python isn't installed. Get it from
  https://www.python.org/downloads/ . **On Windows, during install, tick the box
  that says "Add Python to PATH"** — this is easy to miss and causes exactly the
  "not found" error. Close VS Code and reopen it after installing.

From here on, wherever you see `python`, use whichever word worked for you.

---

## Step 5 — Install the two libraries (one time only)

The project needs two add-on packages, `numpy` and `matplotlib`. Install both
with one command:

```
python -m pip install numpy matplotlib
```

(The `python -m pip` form — rather than just `pip` — guarantees the packages go
to the *same* Python you'll run the scripts with. This avoids the single most
confusing beginner error, where the script later says it can't find numpy.)

You'll see a few lines of download text. When it finishes and you get the prompt
back, you're done. This only has to be done once.

---

## Step 6 — Run the first script

```
python run_bands.py
```

Success looks like exactly one line:

```
Saved figures/bands.png
```

That's it. The script ran, computed the band structures, and saved a picture.

---

## Step 7 — Look at your result

In the left sidebar, open the `figures` folder and click `bands.png`. VS Code
shows images right in the editor. You should see four panels of energy bands.
That's a real result — you just computed and plotted the band structure of a
topological model.

---

## Step 8 — Run the other two

Same pattern:

```
python run_chern.py
```

This prints the Chern numbers (they should be very close to whole numbers like
`+1.0000`) and saves `figures/berry_curvature.png`.

```
python run_phase_diagram.py
```

This saves `figures/phase_diagram.png` — the staircase plot, which is the
headline result of the whole project.

And to confirm everything is correct:

```
python test_models.py
```

Every line should say `PASS` (it checks both the QWZ and the Haldane model).

---

## Step 9 — The fun part: the interactive explorer

```
python explore.py
```

This opens a window with two sliders. Drag the **mass u** slider and watch the
band structure, the Berry curvature map, and the Chern number all update live,
with a banner telling you whether the current phase is topological or trivial.
The **grid N** slider changes the resolution (higher N = smoother, a little
slower).

This opens an interactive desktop window. If you're on a normal laptop it just
works. If no window appears (some remote or WSL-without-display setups lack a
graphical backend), it's not a bug in the project — the static scripts above
still produce the same figures as saved images.

---

## Step 10 — The second model: Haldane

Everything above is the QWZ model (a square lattice). The project also includes
the Haldane model (a honeycomb lattice, the one that started the whole field):

```
python3 run_haldane.py
```

This saves two pictures: the band structure, and the famous 2D phase diagram
showing two topological regions. Then explore it live:

```
python3 explore_haldane.py
```

Three sliders let you change the mass, the phase, and the hopping strength, with
your position tracked on the phase diagram and the Chern number updating live.

---

## Step 11 — The payoff: edge states

```
python3 run_edge_states.py
```

This is the physical point of the whole project. It saves two pictures. The first
shows a topological ribbon whose energy gap is bridged by conducting "edge"
bands, next to a trivial ribbon with a clean, empty gap. The second shows that an
edge state really does hug the edge of the ribbon and fade into the interior. The
script also prints the punchline: the number of edge channels equals the bulk
Chern number.

---

## Quick reference: which "python" word?

| If this works…        | …use this everywhere |
|-----------------------|----------------------|
| `python --version`    | `python`             |
| `python3 --version`   | `python3`            |

Whatever you pick, use the same word for pip and for running scripts:
`python -m pip install …` and `python run_bands.py`, OR
`python3 -m pip install …` and `python3 run_bands.py`.

---

## Troubleshooting

**The terminal is stuck on a line showing `>` or `quote>` or `dquote>`.**You pasted code or an unfinished command into it. Press **Ctrl + C** once to get
back to a normal prompt. Then run files with `python run_bands.py` instead of
pasting their contents.

**`python: command not found` / `Python was not found`.**
Use `python3` instead (see Step 4). If both fail, install Python and tick "Add to
PATH" (Windows), then reopen VS Code.

**`can't open file '...run_bands.py': No such file or directory`.**
The terminal isn't in the project folder. Easiest fix: close the terminal panel,
make sure you did **File > Open Folder** on `chern-calculator` (Step 2), then open
a fresh terminal with **Terminal > New Terminal**. To check where you are, type
`ls` (Mac/Linux) or `dir` (Windows) and confirm you see `run_bands.py` in the
list.

**`ModuleNotFoundError: No module named 'numpy'` (or matplotlib).**
The install went to a different Python than the one running the script. Re-run
the install using the exact same word you run scripts with, e.g.
`python -m pip install numpy matplotlib`. Then try the script again.

**`pip: command not found`.**
Use the `python -m pip …` form instead of bare `pip`. That always works if Python
itself works.

**A window doesn't pop up when I run a script.**
That's expected and fine. These scripts don't open windows — they *save* image
files into the `figures` folder. Open the PNGs from the sidebar to view them.

**Something else / an error I don't recognize.**
Copy the entire terminal output — the command you typed and everything below it —
and paste it in. The exact text is what pinpoints the problem.

---

## Optional: a cleaner setup with a virtual environment

You can skip this entirely; the steps above work fine. A "virtual environment"
just keeps this project's packages separate from the rest of your system, which
is good practice once you do more projects.

```
# create it (one time)
python -m venv .venv
```

Activate it (run this each new terminal session):

- **Mac / Linux:** `source .venv/bin/activate`
- **Windows PowerShell:** `.venv\Scripts\Activate.ps1`
- **Windows Command Prompt:** `.venv\Scripts\activate.bat`

If Windows PowerShell refuses with a message about scripts being disabled, run
this once, then try activating again:

```
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

When the environment is active you'll see `(.venv)` at the start of the prompt.
Then `python -m pip install -r requirements.txt` installs both libraries, and you
run scripts as usual. To leave it, type `deactivate`.
