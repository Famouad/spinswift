/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation
import CGSL

// atomistic code syntax
/*
var a: Atom = Atom(name: PeriodicTable().Z_label("Iron"), type:1, position: Vector3(1,1,1), spin: Vector3(direction:"+y"), g:2)
var b: Atom = Atom(name: PeriodicTable().Z_label("Iron"), type:1, position: Vector3(0,1,2), spin: Vector3(direction:"+y"), g:2)

var atoms: [Atom] = [a,b]

atoms.forEach {
    print(try! $0.jsonify())
    }
*/

// device syntax

/*var MTJ:Stack = Stack()
MTJ.append(Magnetization(name: "Fe", type: 1, position: Vector3(0,0,0), spin: Vector3(direction:"+x"), g:2))
MTJ.append(Magnetization(name: "O", type: 2, position: Vector3(0,0,1), spin: Vector3(0,0,0), g:0))
MTJ.append(Magnetization(name: "Fe", type: 1, position: Vector3(0,0,2), spin: Vector3(direction:"+x"), g:2))

var h: Interaction = Interaction([MTJ[0],MTJ[2]])
h.ZeemanField(Vector3(direction:"+z"), value: 100000)
.UniaxialField(Vector3(direction:"+x"), value: 0.0)
//.ExchangeField(typeI:1,typeJ:1,value:1.0,Rcut:3)
h.Dampening(0.1)

//print(try! h.jsonify())
//let s: Integrate = Integrate(h)
//s.expLs(method:"euler",Δt:0.1)
//print(try! h.atoms[0].jsonify())
//print(Analysis(h.atoms).GetTemperature())
*/

let pulse = LaserExcitation.Pulse(Form: "Gaussian", Fluence: 10.0, Duration: 60E-15, Delay: 0)
let Cp = LaserExcitation.TTM.HeatCapacity(Electron:6E3, Phonon:2.2E6)
let G = LaserExcitation.TTM.Coupling(ElectronPhonon: 2.5E17)
let ttm = LaserExcitation.TTM(EffectiveThickness: 1e-9, InitialTemperature: 300, Damping: 1E-12, HeatCapacity: Cp, Coupling: G)
let laser = LaserExcitation(temperatures: .init(Electron:ttm.InitialTemperature, Phonon:ttm.InitialTemperature), pulse:pulse,ttm:ttm)
print(laser.ComputeInstantPower(time:laser.CurrentTime))

//let timestep = 1E-15
/*for _ in 0...1500 {
    laser.AdvanceTemperaturesGaussian(method:"rk4",Δt:timestep)
    laser.CurrentTime += timestep
    print(laser.CurrentTime,laser.temperatures.Electron,laser.temperatures.Phonon)
<<<<<<< Updated upstream
}
=======
}*/

var data = String()
let mm = Atom.Moments(spin:Vector3(1,0.0,0.0),sigma: Matrix3(1,0,0,0,0,0,0,0,0))
/*
let Fe: Atom = Atom(name: "Fe", type: 1, position: Vector3(0,0,0), ω: -γ.value*1*Vector3(direction:"+z"), moments: mm, g:2)
var time:Double = 0
for _ in 0...400000{
//let s1: Integrate = Integrate(h)
Fe.AdvanceSpindLLB(method: "expIO1", Δt:1e-3, T: 10, α: 0.1)
//print(try! Fe.jsonify())
time+=1e-3
let length=Fe.spin.Norm()
//data+=String(time)+" "+String(Fe.moments.spin.x)+" "+String(Fe.moments.spin.y)+" "+String(Fe.moments.spin.z)+" "+String(length)+"\n"
data+=String(time)+" "+String(Fe.spin.x)+" "+String(Fe.spin.y)+" "+String(Fe.spin.z)+" "+String(length)+"\n"
}
SaveOnFile(data: data, fileName: "dLLBnotemp")
*/

let initials = InitialParam(name:"Ni", type: 1, moments: mm, g: 2.02)
let inisimulation = SimulationProgram.InitializeSimulation(T_initial: 0, T_step: 25, T_final: 250, time_step: 1e-2, stop: 10, α: 0.1)
var Ni: [Atom] = [Atom(),Atom(),Atom(),Atom()]
Ni[0].position = Vector3(0.0, 0.0, 0.0)
Ni[1].position = Vector3(0.0, 0.5, 0.5)
Ni[2].position = Vector3(0.5, 0.0, 0.5)
Ni[3].position = Vector3(0.5, 0.5, 0.0) 
//Define the unit cell atoms (positions in the unit cell). 
let unitCellAtoms: [Atom] = Ni 
// Define the supercell dimensions. 
let supercellDimensions = (x: 1, y: 1, z: 1) 
// Generate the crystal structure. 
let crystalStructure = GenerateCrystalStructure(UCAtoms: unitCellAtoms, supercell: supercellDimensions, LatticeConstant: 0.35, InitParam: initials) 
// Print the positions of atoms in the generated crystal structure. 
//IninitialzeAtomProperties(System: crystalStructure, name: "Ni", type: 1, moments: mm)
/*for atom in crystalStructure {data+=String(atom.name)+" "+String(atom.type)+" "+String(atom.position.x/0.35)+" "+String(atom.position.y/0.35)+" "+String(atom.position.z/0.35)+" "+String(atom.name)+"\n"
print(try! atom.jsonify())}
exit(-1)*/
//let mm = Atom.Moments(spin:Vector3(0,0.435,-0.9),sigma: Matrix3(0,0,0,0,0.19,-0.39,0,-0.39,0.81))
//print(String(crystalStructure.count))
//print(try! Ni.jsonify())
SaveOnFile(data: data, fileName: "struct_NiFcc")
//let distance: Double = ComputeDistance(Boxsize: 0.35*Vector3(5,5,5), PBC: "true", atom1: crystalStructure[0], atom2: crystalStructure[25])
//print(String(distance))
let Boundaries = BoundaryConditions(BoxSize: 0.35*Vector3(1,1,1) ,PBC: "off")

var h: Interaction = Interaction(crystalStructure)
.ExchangeField(typeI:1,typeJ:1,value:13.725/ℏ.value,Rcut:0.25, BCs: Boundaries)
//.ZeemanField(Vector3(direction:"+z"), value: 1)

//print(try! h.jsonify())

let sol: Integrate = Integrate(h)
//sol.Evolve(stop: 1e-1, Δt: 1e-5, method: "rk4", file: "Output_atms", Eqs: "dLLB")

//print(try! h.jsonify())

var p = SimulationProgram(sol)
p.Simulate(Program: "curie_temperature", Initialize: inisimulation)
>>>>>>> Stashed changes
