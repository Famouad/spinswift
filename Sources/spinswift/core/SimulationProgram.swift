/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation
/// This class for managing The simualtion programs of spinswift code
/// All simualtion programs are created and managed in this class, the list of available simulation programs is bellow
/// - Author: Mouad Fattouhi 
/// - Date: 02/09/2024
/// - Version: 0.1

/* List of simulation programs:

   1) Curie_Temperature: This program allows for simulating multispin dynamics with sweept temprerature within a range provided by the user in order to determine the Curie curve. 

   2) Optical_Pulse: This Program simulate the effects of an applied laser pulse on magnetization by conecting a 2 or 3 tempretures model to the
      magnetization dynamics equations.

   3) MacroSpin: This simualtion program allows for simualting a single spin dynamics in presence of user defined magnetic interaction using either sLLG or
      dLLB equations.

   4) ParamgneticSpins: This program allows for simulating a the dynamics of a multispin system in abscence of exhcange interactions. 
*/

class SimulationProgram : Codable {

    var I:Integrate

    init(_ I: Integrate = Integrate()) {
        self.I = I
    }

    struct InitializeSimulation : Codable {
        var T_initial: Double
        var T_step: Double
        var T_final: Double
        var stop: Double
        var time_step: Double
        var α : Double
        
        init(T_initial: Double? = Double(), T_step: Double? = Double(), T_final: Double? = Double(), time_step: Double? = Double(), stop: Double? = Double(), α: Double? = Double()){
            self.T_initial = T_initial!
            self.T_step = T_step! 
            self.T_final = T_final!
            self.time_step = time_step!
            self.stop = stop!
            self.α = α!    
        }
    }

    func Simulate(Program: String? = nil, Initialize: InitializeSimulation) {
        switch Program?.lowercased(){
        case "curie_temperature"? :
        self.CurieTemp(Initialize: Initialize)
        case "optical_pulse"? : break
        //self.OpticalPulse(stop: stop, Δt:Initialize.time_step)
        case "macrospin"? : break
        //self.Macrospin(stop: stop, Δt:Initialize.time_step)
        case "paramagneticspins"? : break 
        //self.Paramagnet(stop: stop, Δt:Initialize.time_step)
        default: print("Program not found. Choose one of the available programs:\ncurie_temperature\noptical_pulse\nmacrospin\nparamagneticspins")
        }
    }

    func CurieTemp(Initialize: InitializeSimulation) {
        var T: Double = Initialize.T_initial
        var content: String = String()
        let stop: Double = Initialize.stop; let Δt : Double = Initialize.time_step; let T_final: Double = Initialize.T_final
        let dT: Double = Initialize.T_step; let α: Double = Initialize.α 
        while (T < T_final) {
            let fnm: String = "Output_CurieTemp_"+String(T)
            I.EvolveRK45(stop: stop, Δt: Δt, fileName: fnm, Equations: "dLLB", Temp: T, alpha: α)
            let m : Vector3 = Analysis(I.h.atoms).GetMagnetization()
            let mnorm : Double = Analysis(I.h.atoms).GetMagnetizationLength()
            content += String(T)+"\t"+String(m.x)+"\t"+String(m.y)+"\t"+String(m.z)+"\t"+String(mnorm)+"\n"
            T+=dT
        } 
        SaveOnFile(data:content, fileName: "Output_CurieTemp_MvsTfile")
    }
}

