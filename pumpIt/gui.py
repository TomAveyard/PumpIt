import PySimpleGUI as sg

# Defines a multline print function and redefines print as this function
# As a result all print statements get printedto the GUI console
def mprint(*args, **kwargs):
    window['-CONSOLE-' + sg.WRITE_ONLY_KEY].print(*args, **kwargs)

print = mprint

# Create custom theme
black_and_white = {'BACKGROUND': '#ffffff',
                'TEXT': '#000000',
                'INPUT': '#ffffff',
                'TEXT_INPUT': '#000000',
                'SCROLL': '#ffffff',
                'BUTTON': ('#000000', '#ffffff'),
                'PROGRESS': ('#ffffff', '#000000'),
                'BORDER': 1,
                'SLIDER_DEPTH': 0,
                'PROGRESS_DEPTH': 0}

# Add to PySimpleGUI themes
sg.theme_add_new('BlackAndWhite', black_and_white)

sg.theme('Black And White')

tab1_layout = [
    [sg.Text("Test 1 ----")],
    [sg.Text("Test 2 -------------")],
    [sg.Text("Test 3")],
    [sg.VPush()]
]

tab2_layout = [
    [sg.Text("Test 1 ----------------")],
    [sg.Text("Test 2 ------")],
    [sg.Text("Test 3 ------------")],
    [sg.Text("Test 4 ------------")],
    [sg.VPush()]
]

input_column = [
    [sg.TabGroup([[sg.Tab("Tab 1", tab1_layout), sg.Tab("Tab 2", tab2_layout)]], expand_x=True, expand_y=True, border_width=0, )]
]

output_column = [
    [sg.Text("Output 1")],
    [sg.Text("Output 2")],
    [sg.Text("Output 3")]
]

layout = [
    [sg.Column(input_column, key="-INPUTCOLUMN-", expand_x=True, expand_y=True, scrollable=True, vertical_scroll_only=True, element_justification="right", ), sg.Column(output_column, key="-OUTPUTCOLUMN-", expand_y=True, expand_x=True, element_justification="c")],
    [sg.MLine(key="-CONSOLE-"+sg.WRITE_ONLY_KEY, size=(40, 10), expand_x=True)],
    [sg.Push(), sg.Button("Calculate", size=(10,1)), sg.Button("Exit")]
]

window = sg.Window("pumpIt Designer", layout, finalize=True, resizable=True)
#window["-INPUTCOLUMN-"].Widget.configure(borderwidth=1, relief=sg.DEFAULT_FRAME_RELIEF)
window["-OUTPUTCOLUMN-"].Widget.configure(borderwidth=1, relief=sg.DEFAULT_FRAME_RELIEF)

print("-"*40)
print("pumpIt Designer v0.1")
print("-"*40)

while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, "Exit"):
        break


window.close()