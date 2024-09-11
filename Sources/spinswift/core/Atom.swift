/*
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/
import Foundation

/// A class for managing atomic properties
/// 
/// An atomic system is a collection of atoms and relationships with them. 
/// - Author: Pascal Thibaudeau
/// - Date: 14/04/2023
/// - Update author: Mouad Fattouhi
/// - Updated: 23/08/2024
/// - Version: 0.1
class Atom : Codable {
    /// Name of the atomic species
    var name : String
    /// Identifier of the atomic type
    var type : Int
    /// Atomic Landé factor in Bohr's magneton unit   
    var g : Double  
    /// Atomic Cartesian position 
    var position : Vector3 
    /// Atomic Cartesian spin direction
    var spin : Vector3 
    /// Atomic pulsation vector around its spin spins
    var ω: Vector3
    /// The first and second moments of dLLB
    var moments: Moments 
    /** The first moment of the dLLB is the averga of the spin <S> as a vector3
        The second moment of the dLLB is the statistical cumulant of the spin <S*S> which is a Matrix3 
    **/
     
    /**
        Initializes a new atom with the provided parts and specifications.

        - Parameters:
            - name: The name of the atom
            - type: The type of the atom
            - position: The Cartesian position of the atom
            - spin: The spin direction of the atom as a unit vector
            - ω: The pulsation vector located on the atom
            - g: The Landé factor of the atom in Bohr's magneton unit
            - moments: The first and second moments of dLLB <S_i> qnd <S_i*S_j> i,j=x,y,z
        
        - Returns: A new atom
    */
    init(name: String? = String(), type: Int? = Int(), position: Vector3? = Vector3(), spin: Vector3? = Vector3(), ω: Vector3? = Vector3(), moments: Moments? = Moments(), g: Double? = Double()){
    //    if name.isEmpty {fatalError("A name has to be provided")}
        self.name = name!
        self.type = type!
        self.position = position!
        self.spin = spin!
        self.ω = ω!
        guard g! >= 0 else {fatalError("g factor must be positive!")}
        self.g = g!
        self.moments=moments!
    }

    struct Moments : Codable {
        var spin : Vector3 
        var sigma: Matrix3

        init(spin: Vector3? = Vector3(), sigma: Matrix3? = Matrix3()) {
            self.spin = spin!
            self.sigma = sigma!
        }

        static func + (a: Moments, b: Moments) -> Moments {
         return Moments(spin: (a.spin+b.spin), sigma: (a.sigma+b.sigma))
        }
        
        static func += (a: inout Moments, b: Moments) {
            var c: Moments = Moments()
            c = a + b
            a = c
        }

        static func * (a: Double, b: Moments) -> Moments {
         return Moments(spin: a*(b.spin), sigma: a*(b.sigma))
        }
    }

    /**
        Compute the atomic spin direction according to a precession motion around its local pulsation vector.

        - Parameters: 
          - method: The type of method used either [euler](euler) or [symplectic](symplectic)
          - Δt: The timestep used (internal unit in ps).

        - Returns: The Cartesian vector pointing the new direction and located on the unit sphere
    */
   func AdvanceSpin(method: String, Δt: Double) {
        var s: Vector3 = Vector3()
        switch method.lowercased() {
        case "euler" :
            s = spin + Δt * (ω × spin)
        case "symplectic" :
            let ω2: Double = ω°ω
            let c: Double = 0.25*Δt*Δt
            let c2: Double = 1.0/(1.0 + c*ω2)
            var s1: Vector3 = c*((2.0*(ω°spin))*ω - ω2*spin)
            s1 += Δt * (ω × spin)
            s = c2 * (spin + s1)
        case "full" :
            let n: Double = ω.Norm()
            let Ω: Vector3 = (1.0/n)*ω
            let ξ: Double = n*Δt
            let χ: Double = Ω°spin
            s = cos(ξ)*spin + sin(ξ)*(Ω × spin) + (χ*(1.0-cos(ξ)))*Ω
        default: break
        }
        spin = s
        spin.Normalize()
    }

    /**
        Compute the atomic spin direction according to a precession motion around its local pulsation vector with both transverse and longitudinal dampings.

        - Parameters: 
          - method: The type of method used either [euler](euler) or [RangeKutta45](rk4)
          - Δt: The timestep used (internal unit in ps).

        - Returns: The Cartesian vector pointing the new direction with a non-conserved norm
    */

    private func RHS(moments: Moments, Temp: Double? = Double(), alpha: Double? = Double()) -> Moments {
        let α: Double = alpha!
        let c: Double = 1/(1+(α*α))
        let β: Double = 1/(k_B.value*Temp!)
        let D: Double = (α)/(β*ℏ.value)
        let spin: Vector3 = moments.spin
        let Σ: Matrix3 = moments.sigma
        let I: Matrix3 = Matrix3(fill:"identity")

        var rhsdLLB: Moments = Atom.Moments()
        
        var A1: Matrix3 = (Σ.Trace()*(ω ⊗ spin)) - (Σ.Transpose()*(ω ⊗ spin))
        A1+=((ω ⊗ spin)*Σ.Transpose()) - ((ω ⊗ spin).Trace()*Σ.Transpose())
        A1+=((ω ⊗ spin) - (ω ⊗ spin).Transpose())*Σ.Transpose()
        A1-=2*(((spin ⊗ spin).Trace()*(ω ⊗ spin)) - ((spin ⊗ spin).Transpose()*(ω ⊗ spin)))  
        let M1: Matrix3 = (ω × Σ.Transpose()) - (α*A1)
        let M2: Matrix3 = (2*Σ.Trace()*I) - (3*(Σ+Σ.Transpose()))

        rhsdLLB.spin=c*((ω × spin) - ((α*Σ.Trace()*ω) - (α*Σ.Transpose()*ω)) - (2*D*c*spin))
        rhsdLLB.sigma=c*(M1 + (D*c*M2) + M1.Transpose())

        return rhsdLLB

    }

    func AdvanceSpindLLB(method: String, Δt: Double, T: Double? = Double(), α: Double? = Double()) {
        switch method.lowercased() {
        case "euler" :
            self.moments += Δt*RHS(moments:self.moments, Temp: T!, alpha: α!)
        case "rk4" :
            let k1: Moments = RHS(moments:self.moments, Temp: T!, alpha: α!)
            let k2: Moments = RHS(moments:self.moments+0.5*Δt*k1, Temp: T!, alpha: α!)
            let k3: Moments = RHS(moments:self.moments+0.5*Δt*k2, Temp: T!, alpha: α!)
            let k4: Moments = RHS(moments:self.moments+Δt*k3, Temp: T!, alpha: α!)
            self.moments += (Δt/6)*(k1+2*k2+2*k3+k4)
        case "expio1" :
            print("Not implemented")
        default: break
        }
        spin = moments.spin
    }
    
    func jsonify() throws -> String {
        let data: Data = try JSONEncoder().encode(self)
        let jsonString: String? = String(data:data, encoding:.utf8) 
        return jsonString!
    } 
}

/*
var A1: Matrix3 = (Σ.Trace()*(ω ⊗ spin)) - (Σ.Transpose()*(ω ⊗ spin))
            A1+=((ω ⊗ spin)*Σ.Transpose()) - ((ω ⊗ spin).Trace()*Σ.Transpose())
            A1+=((ω ⊗ spin) - (ω ⊗ spin).Transpose())*Σ.Transpose()
            A1-=2*(((spin ⊗ spin).Trace()*(ω ⊗ spin)) - ((spin ⊗ spin).Transpose()*(ω ⊗ spin)))  
            var M1: Matrix3 = (ω × Σ.Transpose()) - (α*A1)
            var M2: Matrix3 = (2*Σ.Diag()) - (3*(Σ+Σ.Transpose()))  
            s = spin + Δt*c*((ω × spin) - ((α*Σ.Trace()*ω) - (α*Σ.Transpose()*ω)) - (2*D*c*spin))
            Σ += Δt*c*(M1 + (D*c*M2) + M1.Transpose()) 
*/