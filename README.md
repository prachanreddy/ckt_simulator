Circuit Simulator (Graph Based)

A Python-based symbolic circuit analysis tool that models electrical networks using graph theory and computes branch currents and voltages in the Laplace domain.


* Graph-based modeling using NetworkX
* Cutset and tieset matrix formulation
* Symbolic Laplace-domain analysis (SymPy)
* Supports R, L, C, DC and AC sources
* Unit parsing (k, m, u, etc.)


How to Run


pip install -r requirements.txt<br>
python cktsim.py

Example netlist<br>
AC1 1 0 10<br>
R1 1 2 10<br>
L1 2 3 25m<br>
C1 3 0 100u


Output

* Branch currents I(s)
* Branch voltages V(s)

