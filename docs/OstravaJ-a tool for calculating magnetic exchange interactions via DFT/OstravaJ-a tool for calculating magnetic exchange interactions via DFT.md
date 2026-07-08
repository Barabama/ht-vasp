# OstravaJ: a tool for calculating magnetic exchange interactions via DFT

Jan Priessnitz<sup>a</sup>, Dominik Legut<sup>a</sup>

<sup>a</sup>IT4Innovations, VSB - Technical University of Ostrava, 17. listopadu 2172/15, 708 00 <sup>ˇ</sup> Ostrava-Poruba, Czech Republic

## Abstract

OstravaJ is a Python package for high-throughput calculation of exchange interaction terms in the Heisenberg model for magnetic materials. It uses the total energy diference method, where calculations are based on the total energy of the system in diferent magnetic configurations, calculated by means of density functional theory. OstravaJ can propose a suitable set of magnetic configurations, generate VASP configuration files in cooperation with the user, and read VASP calculation results, which minimizes necessary human interaction. It can also calculate other relevant properties (e. g. MFA and RPA critical temperature, spin-wave stifness) and provide input for various atomistic spin dynamics codes.

We present results for a number of materials from various classes (metals, transition metal oxides), compared to other methods. They show that the total energy diference method is a useful method for exchange interaction calculation from first principles.

Keywords: magnetic materials, exchange interactions, classical Heisenberg model, total energy diference method, density functional theory

## 1. Introduction

Density functional theory (DFT) is a widespread method which can be used to calculate many properties of condensed materials from first principles. Gradual improvements in the accuracy of DFT have recently allowed us to focus on magnetism in materials – usually operating on the meV-eV per atom energy scale.

To better understand magnetism in condensed matter, several microscopic models of magnetism were introduced. A very successful example is the classical Heisenberg model, which reduces the system into a set of rigid spins on atoms influenced by interactions of various origins. This model is simple enough to allow for simulating magnetic systems on greater scales (atomistic spin dynamics, micromagnetics), but still flexible enough to accurately describe the majority of magnetic materials. The most important interaction in the classical Heisenberg model is the exchange interaction, coupling pairs of spins. The strength of the exchange interactions influences many properties of a magnetic material, deeming it suitable or unsuitable for applications. A computational method for calculating exchange interactions from first principles can significantly accelerate the search for new magnetic materials and the development of new electronic devices.

There are several methods for calculating exchange interactions. A popular method is based on the work of Liechtenstein, Katsnelson, Antropov, and Gubanov – the LKAG formula, or the magnetic force theorem [1, 2]. Principally a linear-response method, it calculates exchange interactions by infinitesimally tilting one spin and calculating the resulting energy correction via perturbative theory and Green’s function. A prerequisite for this approach is a DFT calculation of the ground state in a localized basis set. Results in non-local basis sets need to be expressed in terms of Wannier functions [3], for example, via the Wannier90 code [4, 5]. A popular code for automating this calculation is the TB2J package [6]. Time-dependent DFT can also be used [7].

We present the OstravaJ package, which aims to automate an alternative, conceptually simpler method of calculating exchange interaction energies – the total energy diference method [8, 9]. It is based on calculating the total energy of several magnetic configurations (magnetic phases) using DFT and then mapping the total energies onto the Heisenberg Hamiltonian of the system. The result is a system of linear equations with exchange interaction energies as variables. A significant advantage of this method is that it does not assume any requirements on the basis set of the DFT calculation. On the other hand, some problems may arise when calculating the total energy of the excited magnetic configurations. Furthermore, finding suitable magnetic configurations to calculate long-range interaction energies becomes increasingly dificult. OstravaJ introduces a novel scheme for suitable magnetic configuration selection which aims to tackle this challenge.

## 2. Methods and implementation

## 2.1. Total energy diference method

Let us have an infinite periodic system of magnetic ions with atomic magnetic moments pointing in the direction defined by the spin vectors $\vec { S } _ { i }$ . The system can be described by the following classical Heisenberg Hamiltonian:

$$
\mathcal {H} = \mathcal {H} _ {0} - \sum_ {<   i j >} J _ {i j} ^ {\prime} \vec {S} _ {i} \cdot \vec {S} _ {j} \qquad \vec {S} _ {i} \in \mathbb {R} ^ {3}, | \vec {S} _ {i} | = 1
$$

where $\mathcal { H } _ { 0 }$ denotes the non-magnetic part of the energy and $J _ { i j } ^ { \prime }$ is the strength (energy) of the exchange interaction between spins $\vec { S } _ { i }$ and $\vec { S } _ { j }$ . The summation goes over all spin pairs once.

Due to the spatial symmetry of the system, we can assume that certain exchange interactions must be equal. We can thus sort all spin pairs $( i , j )$ equivalence classes $X _ { k }$ based on the exchange interaction energies $J _ { i j } ^ { \prime }$ and define $J _ { k } \mathrm { ~ - ~ a ~ }$ representative of k-th exchange interaction equivalence class, such that $\forall ( i , j ) \in X _ { k } : J _ { i j } ^ { \prime } = J _ { k }$

The Heisenberg Hamiltonian of the system can then be expressed as

$$
\mathcal {H} _ {\mathrm{Heis}} = \sum_ {k} J _ {k} c _ {k}
$$

$$
c _ {k} = \sum_ {<   i j > \in X _ {k}} \vec {S} _ {i} \cdot \vec {S} _ {j}
$$

In an infinite system, the number of equivalence classes is also infinite. OstravaJ introduces a distance cutof $d _ { \mathrm { m a x } }$ , where only spin pairs separated by distance smaller than $d _ { \mathrm { m a x } }$ are considered. This limits the number of equivalence classes to $N _ { k }$ and all other exchange interaction energies are considered zero. The Heisenberg Hamiltonian then becomes

$$
\mathcal {H} _ {\mathrm{Heis}} = \sum_ {k <   N _ {k}} J _ {k} c _ {k}
$$

Let us say that we have defined $N _ { l }$ concrete magnetic configurations and calculated their total energies $E ^ { ( l ) }$ . Let us denote the final (converged) direction of i-th spin in l-th magnetic configuration as $\vec { S } _ { i } ^ { ( l ) }$ and k-th coeficient as $c _ { k } ^ { ( l ) }$

$$
c _ {k} ^ {(l)} = \sum_ {<   i j > \in X _ {k}} \vec {S} _ {i} ^ {(l)} \cdot \vec {S} _ {j} ^ {(l)}
$$

For each magnetic configuration, we can construct an expression for the energy $E ^ { ( l ) }$ with $J _ { k }$ as free variables

$$
E ^ {(l)} = E _ {0} + \sum_ {k <   N _ {k}} J _ {k} c _ {k} ^ {(l)}
$$

where $E _ { 0 }$ is the non-magnetic part of the total energy.

We can now assemble a system of linear equations with $N _ { l }$ equations and $N _ { k } + 1$ free variables $( J _ { 1 } , . . . , J _ { N _ { k } }$ and $E _ { 0 } )$ which can be represented by the following matrix for the case of $N _ { k } = 3$ and $N _ { l } = 4$

$$
\left( \begin{array}{c c c c c} c _ {1} ^ {(1)} & c _ {2} ^ {(1)} & c _ {3} ^ {(1)} & 1 & E ^ {(1)} \\ c _ {1} ^ {(2)} & c _ {2} ^ {(2)} & c _ {3} ^ {(2)} & 1 & E ^ {(2)} \\ c _ {1} ^ {(3)} & c _ {2} ^ {(3)} & c _ {3} ^ {(3)} & 1 & E ^ {(3)} \\ c _ {1} ^ {(4)} & c _ {2} ^ {(4)} & c _ {3} ^ {(4)} & 1 & E ^ {(4)} \end{array} \right)\tag{1}
$$

Solving this matrix yields a vector of values $\left( J _ { 1 } , J _ { 2 } , J _ { 3 } , E _ { 0 } \right)$ . The system of linear equations can then be solved exactly if $N _ { l } = N _ { k } + 1$ or by the least-squares method if $N _ { l } > N _ { k } + 1$ . Note that a unique solution exists only if the matrix has full rank and the vectors represented by the matrix rows form a basis spanning the full $( N _ { k } + 1 )$ )-dimensional space. In the opposite case, a non-unique solution for $J _ { k }$ carries little physical information. This requirement can be satisfied by choosing a suitable set of magnetic exchange interactions, which is described later in this article.

The solution of this system consists of the exchange interactions $J _ { 1 } , . . . , J _ { N _ { k } }$ and the non-magnetic part of total energy $E _ { 0 }$ . When $N _ { l } > N _ { k } + 1$ , the error of the least-squares fitting is also available.

The non-magnetic total energy $E _ { 0 }$ can be further used to easily calculate the critical temperature of the system in the mean-field approximation (MFA). Let $E _ { \mathrm { g s } }$ be the total energy of the ground state and $N _ { s }$ be the number of magnetic sites in the unit cell. The MFA critical temperature can then be calculated as:

$$
T _ {C} ^ {M F A} = \frac {\Delta E}{3 k _ {b} N _ {s}} = \frac {E _ {0} - E _ {g s}}{3 k _ {b} N _ {s}}\tag{2}
$$

## 2.2. Choosing a set of magnetic configurations

First task within the total energy diference method is to find a suitable set of magnetic configurations which leads to a system of linear equations with a unique solution, as described in the section above. This task depends on the interaction distance cutof and, in turn, the number of unique exchange interactions to take into account (denoted $N _ { k } ) \mathrm { ~ - ~ a ~ }$ parameter configured by the user beforehand. The configuration set must include at least $N _ { k } + 1$ configurations which yield linearly independent Heisenberg Hamiltonians.

Searching for magnetic configurations within the primitive unit cell of the material in order to find the suitable set is usually unsuccessful. In the extreme case of a one-atom unit cell, only the ferromagnetic configuration is accessible, and the configuration set can never be assembled. Before the search, the unit cell must be extended to a supercell containing $N _ { x } \times N _ { y } \times$ $N _ { z }$ unit cells. The supercell must be large enough so that a suitable set of $N _ { k } + 1$ independent configurations exists, but not too large so that the memory and computational requirements of the subsequent calculations are reasonable. As of now, the suitable interaction cutof and supercell size need to be configured by the user. OstravaJ will be able to select the best supercell size automatically in the future.

For a unit cell with $N _ { s }$ magnetic sites, the number of possible collinear magnetic configurations $N _ { c }$ grows exponentially as $N _ { c } = 2 ^ { N _ { s } }$ . This makes brute-force traversal of the whole configuration space very impractical. The novelty of OstravaJ lies in the special algorithm for generating the suitable configuration set, able to search through space consisting of hundreds of magnetic sites and more than $2 ^ { 1 0 0 }$ individual configurations.

OstravaJ builds the configuration set incrementally, adding configurations into the set one by one. A configuration can be added only if it fulfills the following criteria:

1. new configuration is linearly-independent to the current set of configurations – if it were added, the rank of the resulting system of equations matrix would increase by 1.

2. equivalent magnetic sites must have identical local magnetic surrounding (must be chemically equivalent)

Let us look more closely at the second requirement of the identical local magnetic surrounding. This can be expressed in terms of the Heisenberg Hamiltonian. First, let us divide the Hamiltonian into contributions made by individual magnetic sites, denoted as $H _ { i }$

$$
\mathcal {H} = \sum_ {<   i j >} J _ {i j} \vec {S} _ {i} \cdot \vec {S} _ {j} = \frac {1}{2} \sum_ {i} H _ {i}\tag{3}
$$

$$
H _ {i} = \sum_ {j} J _ {i j} \vec {S _ {i}} \cdot \vec {S _ {j}} = \vec {S _ {i}} (\sum_ {j} J _ {i j} \vec {S _ {j}})\tag{4}
$$

The second requirement says that individual contributions $H _ { i }$ of all magnetic sites of equal type must be equal

$$
\forall i, j: H _ {i} = H _ {j}\tag{5}
$$

![](images/80bbc88cd31a22c0993b5781da1dd3d86e5068815ef132ddc88172e73e0d79ef.jpg)  
Figure 1: Example of the identical magnetic surrounding criterion. Blue circles show magnetic sites with spin up, red circles show spin down. On the left, the criterion is fulfilled. On the right, the criterion is broken – various sites have various surroundings.

Figure 1 shows an example of the identical magnetic surrounding criterion on a 2D square lattice. The configuration on the left fulfills the criterion – all sites have 4 nearest neighbors pointing antiparallel, 4 next-nearest neighbors pointing parallel, etc. The configuration on the right does not fulfill the criterion – there are sites with 1 nearest neighbor parallel and 3 antiparallel, sites with 2 parallel and 2 antiparallel, and sites with 3 parallel and 1 antiparallel. Various magnetic sites have diferent local magnetic surroundings.

The performance in searching through the configuration space achieved by OstravaJ lies in the fact that it traverses the space recursively and checks whether the identical magnetic surrounding criterion is reachable even in incomplete configurations. In case it is not reachable, it aborts that branch of the traversal early and saves a significant amount of steps. An example of recursive traversal and early abortion is shown in Figure 2.

![](images/5472c8ec05bbcabae78e4506de6e0a14d96df525f3f98a23e7966acebed2c3ff.jpg)  
Figure 2: Example of the recursive traversal and early abortion. In the 7th step, the algorithm is already certain that there is no way to fill the rest of the spins and fulfill the identical magnetic surrounding criterion and it aborts that branch of the search.

## 2.3. Installation

OstravaJ is available at https://code.it4i.cz/jpriessnitz/ostravaj as a Git repository. The package itself does not require explicit installation, but it depends on several Python packages, as listed in the requirements.txt file. Dependencies can be installed via

pip3 install -r requirements.txt

OstravaJ can then be used through

<path-to-OstravaJ>/OstravaJ.sh <subcommand> <args>

Apart from OstravaJ, users must also install a DFT software (currently, only VASP can be used) to conduct complete calculations. However, OstravaJ does not invoke the DFT software directly, so it does not have to be configured for cooperation.

Further instructions on how to use OstravaJ are available in the README.md file in the source code.

## 2.4. Workflow overview

The OstravaJ tool is a Python package which aims to automate the total energy diference method for calculating exchange interactions. It requires external DFT software to conduct the total energy calculations. Currently, the only supported DFT code is VASP [10, 11], but we plan to implement support for other popular DFT codes in the future.

Let us now look at the whole calculation workflow from the user perspective. We will look at each step in more detail in later sections of this article.

Before the calculation begins, user must provide several things:

• crystal structure in the form of a POSCAR file and information about which ions are to be considered magnetic

• size of the magnetic unit cell

• distance cutof or a selection of exchange interactions

• base spin

• configuration files for the VASP DFT calculation (INCAR, KPOINTS, POTCAR, ...)

OstravaJ then begins to search for a suitable set of magnetic configurations based on the crystal structure, magnetic unit cell size and distance cutof. Either it succeeds with the search and outputs a number of VASP configurations, or it fails to find a suitable set which would lead to solving a system with the desired set of exchange interactions. The number of suitable configurations depends on the magnetic unit cell size – more complicated magnetic configurations are accessible in larger magnetic unit cells. In general, the larger the number of exchange interactions, the larger the magnetic unit cell must be.

As a next step, the user must launch the VASP calculations generated by OstravaJ.

After converging all DFT calculations, OstravaJ reads the converged magnetic moments and total energies, assembles the system of linear equations and solves for exchange interaction energies. It also calculates some experimentally important quantities such as critical temperature in the mean-field approximation (MFA) and random-phase approximation (RPA) and spinwave stifness.

Optionally, OstravaJ can also prepare system configuration for atomistic spin dynamics simulations and calculate the magnetization vs. temperature curve. It automatically determines the spin-dynamics critical temperature which is generally considered more precise compared to MFA and RPA (usually $T _ { \mathrm { C , R P A } } < T _ { \mathrm { C , A S D } } < T _ { \mathrm { C , M F A } }$ . Currently, only UppASD code [12] is supported for spin dynamics calculations.

## 3. Results

We showcase the capabilities of OstravaJ on calculations of several simple structures. The results are compared to both computational and experimental data from other sources to see the accuracy of the total energy diference method. All DFT calculations were conducted with the Vienna Ab initio Simulation Package (VASP) software (version 6.4). VASP is a plane-wave basis set implementation [10, 11] of the Density Functional Theory within the Projector Augmented Wave (PAW) approximation [13].

## 3.1. Transition metal oxides

Firstly, we calculate magnetic exchange interactions in transition metal oxides CoO, NiO and MnO. These materials all have rock-salt structure which has face-centered cubic lattice and a unit cell containing one atom of oxygen and one atom of transition metal. They are dielectric materials, meaning that long-range exchange interactions are negligible.

As a first step, we provide the material structure files to OstravaJ and set the exchange interaction cutof to 6 nearest exchange interactions, meaning that OstravaJ needs to find at least 7 unique (linearly independent) magnetic configurations. The supercell size was set to $4 \times 4 \times 2$ unit cells – OstravaJ did not find enough independent configurations in smaller supercells.

Total energies of all configurations were then calculated using VASP with calculation parameters shown in table Table 1. We used the generalizedgradient-approximation with the Perdew-Burke-Ernzerhof (PBE) exchangecorrelation functional [14] with phenomenological Hubbard U via the Dudarev approach [15] to capture the electron correlations. In addition, the nearest two exchange interactions in MnO and NiO were recalculated with the HSE06 hybrid exchange-correlation functionals [16].

<table><tr><td></td><td>MnO(PBE+U)</td><td>MnO(HSE06)</td><td>NiO(PBE+U)</td><td>NiO(HSE06)</td><td>CoO(PBE+U)</td></tr><tr><td>lattice parameter [Å]</td><td>4.45</td><td>4.45</td><td>4.20</td><td>4.20</td><td>4.25</td></tr><tr><td>energy cutoff [eV]</td><td>700</td><td>500</td><td>500</td><td>500</td><td>700</td></tr><tr><td>Hubbard U-J [eV]</td><td>5.0</td><td>-</td><td>5.0</td><td>-</td><td>4.0</td></tr><tr><td>supercell</td><td>4 × 4 × 2</td><td>2 × 2 × 2</td><td>4 × 4 × 2</td><td>2 × 2 × 2</td><td>4 × 4 × 2</td></tr><tr><td>k-grid</td><td>2 × 2 × 4</td><td>6 × 6 × 6</td><td>2 × 2 × 4</td><td>6 × 6 × 6</td><td>2 × 2 × 4</td></tr><tr><td>energy convergence criterion [eV]</td><td> $10^{-7}$ </td><td> $10^{-7}$ </td><td> $10^{-7}$ </td><td> $10^{-7}$ </td><td> $10^{-7}$ </td></tr><tr><td>no. of valence electrons (M+O)</td><td>7+6</td><td>13+6</td><td>16+6</td><td>16+6</td><td>17+6</td></tr></table>

Table 1: Important parameters for the transition-metal-oxide DFT calculations.

![](images/9774d6583efa33c9e817972c2f8dcd7ae582a3eeeca96294134a4f20cd2de2ac.jpg)  
Figure 3: Magnetic exchange interactions of a) MnO: calculated with PBE+U and HSE06 functionals, from inelastic neutron scattering [17], from thermodynamic data [18]; b) NiO: calculated with PBE+U and HSE06 functionals, from inelastic neutron scattering [19], from thermodynamic data [20]; and CoO: calculated with PBE+U approach, calculated with $\mathrm { L D A + U }$ in [21].

Figure 3 shows magnetic exchange interactions calculated with OstravaJ, compared to other computational or experimental works. Unlike in other works, we were able to calculate not only the first 2 nearest exchange interactions, but also 4 further longer-range interactions, and confirm that they are indeed very small. We can see a larger discrepancy in MnO between our calculation employing the $\mathrm { P B E + U }$ approach and the experimental figures.

However, this can be attributed to an inaccuracy in the PBE functional. The match to experiment was significantly improved after using the HSE06 hybrid functional, which is much more accurate, albeit very computationally expensive.

<table><tr><td></td><td>MnO</td><td>NiO</td><td>CoO</td></tr><tr><td>OstravaJ - PBE+U - MFA</td><td>83 K</td><td>444 K</td><td>268 K</td></tr><tr><td>expt. Js + MC</td><td>85 K [17]</td><td>340 K [19]</td><td>210 K [21]</td></tr></table>

Table 2: Calculated N´eel temperatures for MnO, NiO and CoO using a) exchange interactions from this work and mean-field approximation, b) experimental exchange interactions with Monte Carlo calculations

Table 2 shows N´eel temperatures calculated from the exchange interactions using mean-field theory, compared to other sources. For NiO and CoO, the mean-field temperature is significantly larger than the Monte Carlo one. This is most likely due to a well-known overestimation of critical temperatures by the mean-field approximation method.

## 3.2. Metals

The capability of OstravaJ to calculate long-range exchange interactions can be very well shown on a few important metallic materials, which do not have an insulating gap and, unlike in TMOs, the longer-range interactions are non-negligible. In this section, we calculate exchange interactions in hexagonal-closely-packed Cobalt (hcp Co), face-centered-cubic Cobalt (fcc Co), face-centered-cubic Nickel (fcc Ni) and body-centered-cubic Iron (bcc Fe). These materials have been chosen for the availability of published data for comparison.

For each system, we manually selected the exchange interaction cutof and the size of the supercell such that OstravaJ was able to find a suficient amount of independent configurations and the supercell was not too large, so that the calculations were not too computationally complex. The selection is shown in Table 3

<table><tr><td></td><td>hcp Co</td><td>fcc Co</td><td>fcc Ni</td><td>bcc Fe</td></tr><tr><td>no. of nearest exchange interactions</td><td>8</td><td>6</td><td>6</td><td>6</td></tr><tr><td>distance cutoff [Å]</td><td>5.9</td><td>6.1</td><td>6.1</td><td>5.7</td></tr><tr><td>supercell size</td><td>4 × 2 × 2</td><td>4 × 4 × 2</td><td>4 × 4 × 2</td><td>4 × 4 × 2</td></tr></table>

Table 3: Choice of exchange interaction cutof and supercell size for calculations of metals.

Total energies of all magnetic configurations were then calculated by VASP with the PBE exchange-correlation functional and parameters shown in Table 4.

<table><tr><td></td><td>hcp Co</td><td>fcc Co</td><td>fcc Ni</td><td>bcc Fe</td></tr><tr><td>lattice parameter [Å]</td><td>2.48</td><td>3.51</td><td>4.20</td><td>4.20</td></tr><tr><td>energy cutoff [eV]</td><td>400</td><td>500</td><td>400</td><td>400</td></tr><tr><td>k-grid</td><td> $8 \times 15 \times 8$ </td><td> $10 \times 10 \times 20$ </td><td> $8 \times 8 \times 16$ </td><td> $6 \times 6 \times 12$ </td></tr><tr><td>energy convergence criterion [eV]</td><td> $10^{-7}$ </td><td> $10^{-6}$ </td><td> $10^{-6}$ </td><td> $10^{-7}$ </td></tr><tr><td>no. of valence electrons</td><td>9</td><td>17</td><td>10</td><td>8</td></tr></table>

Table 4: Important parameters for the DFT calculations of metals.

The main results of the calculations, the exchange interaction energies, are shown in Figure 4 for hcp Co and fcc Co and Figure 5 for fcc Ni and bcc Fe. Energies calculated via the TB-LMTO (tight-binding linear-mufin-tinorbital) method [22, 2] are shown for comparison.

![](images/5962332843321317402edb3eb3e82209fc1def5fcc9b29b2c72802ec5c11687f.jpg)  
Figure 4: Exchange interactions for a) hcp Co: calculated using OstravaJ and PBE+U approach, via TB-LMTO method [22], b) fcc Co: calculated using OstravaJ and PBE+U approach, via TB-LMTO method [2].

The results from this work and from the TB-LMTO method for the Cobalt systems are qualitatively in agreement – the sign and the magnitude of the nearest interaction energies are comparable, and the energies diminish with distance in a similar fashion.

![](images/69f342bf058b72053d94e7a4c3719beaa621ec08eb1eadce16e852a668d36f2d.jpg)

![](images/71ff48a446dc8fc0772cd05275316c17add6f8111dd5d0e5e239e6dac23ad357.jpg)  
Figure 5: Exchange interactions for a) fcc Ni: calculated using OstravaJ and PBE+U approach, via TB-LMTO method [2], b) bcc Fe: calculated using OstravaJ and PBE+U approach, via TB-LMTO method [2].

Results for fcc Ni and bcc Fe, unfortunately, show a significant deviation from the TB-LMTO method.

<table><tr><td colspan="2"></td><td> $T_{C}$  [K]</td><td> $D$  [meV · Å2]</td></tr><tr><td>hcp Co</td><td>PBE+UTB-LMTO [2]</td><td>1442 (MFA), 1281 (RPA)1673 (MFA)</td><td>597-</td></tr><tr><td>fcc Co</td><td>PBE+UTB-LMTO [2]exp.</td><td>1442 (MFA), 1281 (RPA)1645 (MFA), 1311 (RPA)1388 - 1398</td><td>597663580 [23], 510 [24]</td></tr><tr><td>fcc Ni</td><td>PBE+UTB-LMTO [2]exp.</td><td>416 (MFA), 367 (RPA)397 (MFA), 350 (RPA)624-631</td><td>597756555 [25], 422 [23]</td></tr><tr><td>bcc Fe</td><td>PBE+UTB-LMTO [2]exp.</td><td>1530 (MFA), 708 (RPA)1414 (MFA), 950 (RPA)1044-1045</td><td>103250280 [23], 330 [24]</td></tr></table>

Table 5: Comparison of critical temperature and spin-wave stifness between this work, TB-LMTO method [2] and experimental results.

Exchange interaction energies are by themselves dificult to measure experimentally. However, they can be used to calculate critical temperature and spin-wave stifness $( T _ { \mathrm { C } } , D )$ which are experimentally accessible quantities. Table 5 shows the comparison of $T _ { \mathrm { C } }$ and D among values calculated from OstravaJ method, TB-LMTO method, and measured experimentally. The majority of the figures agree within tens of percent – an adequate accuracy in the context of ab-initio methods. One exception is the spin-wave stifness in bcc Fe – the OstravaJ method reports roughly $D _ { \mathrm { O J } } \approx 1 0 0 \ \mathrm { m e V \cdot \mathring { A } ^ { 2 } }$ , while the experimental value is roughly $D _ { \mathrm { e x } } \approx 3 0 0 \ \mathrm { m e V \cdot \mathring { A } ^ { 2 } }$

## 4. Conclusion

The total energy diference method is a widespread method for calculating exchange interaction energies for a wide variety of magnetic materials. However, a challenging part of this method is the search for and selection of suitable magnetic configurations. We propose a novel method of traversing the magnetic configuration space and overcoming this challenge.

We introduce the OstravaJ code, which automates most parts of the calculation via the total energy diference method, including the search for suitable magnetic configurations. It can thus be easily used in the field of highthroughput computational materials science.

Furthermore, the total energy diference method and the OstravaJ code were benchmarked on a set of several sample materials. We have shown that the total energy diference method is a viable computational method and gives reasonable results. However, the results are not in full agreement with calculations via the TB-LMTO method, which uses the LKAG formula. This can be attributed to the fact that this method is based on perturbing the magnetic ground state, while the total energy diference method also takes into account excited magnetic configurations. Thus, for materials which do not perfectly fit the Heisenberg model picture, the results are going to difer.

## 5. Acknowledgments

The authors acknowledge grant No. 22-35410K by Czech Science Foundation. This work was supported by the Ministry of Education, Youth and Sports of the Czech Republic through the e-INFRA CZ (ID:90254).

## References

[1] A. Liechtenstein, M. Katsnelson, V. Antropov, V. Gubanov, Local spin density functional approach to the theory of exchange interactions in ferromagnetic metals and alloys, Journal of Magnetism and Magnetic Materials 67 (1) (1987) 65–74. doi:https://doi.org/10.1016/0304-8853(87)90721-9. URL https://www.sciencedirect.com/science/article/pii/ 0304885387907219

[2] M. Pajda, J. Kudrnovsk´y, I. Turek, V. Drchal, P. Bruno, Ab initio calculations of exchange interactions, spin-wave stifness constants, and Curie temperatures of Fe, Co, and Ni, Physical Review B 64 (17) (2001) 174402. doi:10.1103/PhysRevB.64.174402. URL https://link.aps.org/doi/10.1103/PhysRevB.64.174402

[3] D. M. Korotin, V. V. Mazurenko, V. I. Anisimov, S. V. Streltsov, Calculation of exchange constants of the heisenberg model in plane-wave-based methods using the green’s function approach, Phys. Rev. B 91 (2015) 224405. doi:10.1103/PhysRevB.91.224405. URL https://link.aps.org/doi/10.1103/PhysRevB.91.224405

[4] G. Pizzi, V. Vitale, R. Arita, S. Bl¨ugel, F. Freimuth, G. G´eranton, M. Gibertini, D. Gresch, C. Johnson, T. Koretsune, J. Iba˜nez-Azpiroz, H. Lee, J.-M. Lihm, D. Marchand, A. Marrazzo, Y. Mokrousov, J. I. Mustafa, Y. Nohara, Y. Nomura, L. Paulatto, S. Ponc´e, T. Ponweiser, J. Qiao, F. Th¨ole, S. S. Tsirkin, M. Wierzbowska, N. Marzari, D. Vanderbilt, I. Souza, A. A. Mostofi, J. R. Yates, Wannier90 as a community code: new features and applications, Journal of Physics: Condensed Matter 32 (16) (2020) 165902. doi:10.1088/1361-648X/ab51ff. URL https://dx.doi.org/10.1088/1361-648X/ab51ff

[5] A. A. Mostofi, J. R. Yates, Y.-S. Lee, I. Souza, D. Vanderbilt, N. Marzari, wannier90: A tool for obtaining maximally-localised wannier functions, Computer Physics Communications 178 (9) (2008) 685–699. doi:https://doi.org/10.1016/j.cpc.2007.11.016. URL https://www.sciencedirect.com/science/article/pii/ S0010465507004936

[6] X. He, N. Helbig, M. J. Verstraete, E. Bousquet, Tb2j: A python package for computing magnetic interaction parameters, Computer Physics Communications 264 (2021) 107938. doi:https://doi.org/10.1016/j.cpc.2021.107938. URL https://www.sciencedirect.com/science/article/pii/ S0010465521000679

[7] T. Gorni, O. Baseggio, P. Delugas, S. Baroni, I. Timrov, turbomagnon – a code for the simulation of spin-wave spectra using the liouville-lanczos approach to time-dependent density-functional perturbation theory, Computer Physics Communications 280 (2022) 108500. doi:https://doi.org/10.1016/j.cpc.2022.108500. URL https://www.sciencedirect.com/science/article/pii/ S0010465522002193

[8] T. Steenbock, C. Herrmann, Toward an automated analysis of exchange pathways in spin-coupled systems, Journal of Computational Chemistry 39 (2) (2018) 81–92. arXiv:https: //onlinelibrary.wiley.com/doi/pdf/10.1002/jcc.25081, doi:https://doi.org/10.1002/jcc.25081. URL https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc. 25081

[9] H. Xiang, C. Lee, H.-J. Koo, X. Gong, M.-H. Whangbo, Magnetic properties and energy-mapping analysis, Dalton Transactions 42 (4) (2012) 823–853, publisher: The Royal Society of Chemistry. doi:10.1039/C2DT31662E. URL https://pubs.rsc.org/en/content/articlelanding/2013/ dt/c2dt31662e

[10] G. Kresse, J. Furthm¨uller, Eficient iterative schemes for ab initio totalenergy calculations using a plane-wave basis set, Phys. Rev. B 54 (1996) 11169–11186. doi:10.1103/PhysRevB.54.11169. URL https://link.aps.org/doi/10.1103/PhysRevB.54.11169

[11] G. Kresse, D. Joubert, From ultrasoft pseudopotentials to the projector augmented-wave method, Phys. Rev. B 59 (1999) 1758–1775. doi:10. 1103/PhysRevB.59.1758. URL https://link.aps.org/doi/10.1103/PhysRevB.59.1758

[12] B. Skubic, J. Hellsvik, L. Nordstr¨om, O. Eriksson, A method for atomistic spin dynamics simulations: implementation and examples, Journal of Physics: Condensed Matter 20 (31) (2008) 315203. doi: 10.1088/0953-8984/20/31/315203. URL https://dx.doi.org/10.1088/0953-8984/20/31/315203

[13] P. E. Bl¨ochl, Projector augmented-wave method, Phys. Rev. B 50 (1994) 17953–17979. doi:10.1103/PhysRevB.50.17953. URL https://link.aps.org/doi/10.1103/PhysRevB.50.17953

[14] J. P. Perdew, K. Burke, M. Ernzerhof, Generalized gradient approximation made simple, Phys. Rev. Lett. 77 (1996) 3865–3868. doi: 10.1103/PhysRevLett.77.3865. URL https://link.aps.org/doi/10.1103/PhysRevLett.77.3865

[15] S. L. Dudarev, G. A. Botton, S. Y. Savrasov, C. J. Humphreys, A. P. Sutton, Electron-energy-loss spectra and the structural stability of nickel oxide: An lsda+u study, Phys. Rev. B 57 (1998) 1505–1509. doi: 10.1103/PhysRevB.57.1505. URL https://link.aps.org/doi/10.1103/PhysRevB.57.1505

[16] J. Heyd, G. E. Scuseria, M. Ernzerhof, Hybrid functionals based on a screened Coulomb potential, The Journal of Chemical Physics 118 (18)

(2003) 8207–8215, eprint: https://pubs.aip.org/aip/jcp/articlepdf/118/18/8207/19093575/8207 1 online.pdf. doi:10.1063/1. 1564060. URL https://doi.org/10.1063/1.1564060

[17] M. Kohgi, Y. Ishikawa, Y. Endoh, Inelastic neutron scattering study of spin waves in mno, Solid State Communications 11 (2) (1972) 391–394. doi:https://doi.org/10.1016/0038-1098(72)90255-4. URL https://www.sciencedirect.com/science/article/pii/ 0038109872902554

[18] M. E. Lines, E. D. Jones, Antiferromagnetism in the face-centered cubic lattice. ii. magnetic properties of mno, Phys. Rev. 139 (1965) A1313– A1327. doi:10.1103/PhysRev.139.A1313. URL https://link.aps.org/doi/10.1103/PhysRev.139.A1313

[19] M. T. Hutchings, E. J. Samuelsen, Measurement of spin-wave dispersion in nio by inelastic neutron scattering and its relation to magnetic properties, Phys. Rev. B 6 (1972) 3447–3461. doi:10.1103/PhysRevB.6.3447. URL https://link.aps.org/doi/10.1103/PhysRevB.6.3447

[20] R. Shanker, R. A. Singh, Analysis of the exchange parameters and magnetic properties of nio, Phys. Rev. B 7 (1973) 5000–5005. doi: 10.1103/PhysRevB.7.5000. URL https://link.aps.org/doi/10.1103/PhysRevB.7.5000

[21] T. Archer, R. Hanafin, S. Sanvito, Magnetism of coo polymorphs: Density functional theory and monte carlo simulations, Phys. Rev. B 78 (2008) 014431. doi:10.1103/PhysRevB.78.014431. URL https://link.aps.org/doi/10.1103/PhysRevB.78.014431

[22] I. Turek, J. Kudrnovsk´y, V. Drchal, P. Bruno, S. Bl¨ugel, Ab initio theory of exchange interactions in itinerant magnets, physica status solidi (b) 236 (2) (2003) 318–324. doi:10.1002/pssb.200301671. URL https://onlinelibrary.wiley.com/doi/10.1002/pssb. 200301671

[23] R. Pauthenet, Experimental verification of spin-wave theory in high fields (invited), Journal of Applied Physics 53 (11) (1982) 8187–8192.

doi:10.1063/1.330287. URL https://doi.org/10.1063/1.330287

[24] G. Shirane, V. J. Minkiewicz, R. Nathans, Spin Waves in 3d Metals, Journal of Applied Physics 39 (2) (1968) 383–390. doi:10.1063/1. 2163453. URL https://doi.org/10.1063/1.2163453

[25] H. A. Mook, J. W. Lynn, R. M. Nicklow, Temperature Dependence of the Magnetic Excitations in Nickel, Physical Review Letters 30 (12) (1973) 556–559, publisher: American Physical Society. doi:10.1103/ PhysRevLett.30.556. URL https://link.aps.org/doi/10.1103/PhysRevLett.30.556