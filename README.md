# Nuclear Physics Excitation Level Scheme Drawer

#### This is my very first project. It is an elementary Python program intended to help physicists draw nuclear excitation-level schemes.
\
\
To call the program move to **src** dir and use **$ python3 draw.py**. This will open a GUI that allows the user to upload the **.csv** files containing the level scheme info and to set everything before drawing.
If the upload went well, then the first two fields will become green, otherwise, they will be highlighted in red and a popup window will tell you which exception occurred.


The program requires two inputs **.csv** files:

- The first one is the **transitions** file that needs to be filled as follow *(Transition_energy, Starting_level, Ending_level, Spin_parity, Color)*
- The second one is the **levels** file that needs to be filled as follow *(Level_energy, Spin_Parity, Level_Position , Color)*
    
N.B. 
- When saving, the currently drawn figure will be saved with the graphical ratio that are currently used to display it. If you want to resize the figure you can just resize the window.
- You can draw 'White' levels at the top of the level scheme to make different level schemes with the same energy-to-size ratio and make them comparable.

GUI preview
![GUI preview](https://github.com/MassiGitRep/Level-Scheme-Drawer-for-Nuclear-Physics/blob/main/images/PyQT5_GUI.png)

\
\
Any suggestion is warmly welcomed!

Cheers,
\
Max
