import math

class Wing:
    def init(self, wing_index, wing_length, mean_chord, r33, r22, r11, r00, z_cp2, z_cp1, z_cp0, z_rd, shoulder_width, stroke_plane_offset):
        self.sign_ = pow(-1,wing_index);
        self.R_w_ = wing_length
        self.d_0_ = shoulder_width
        self.d_s_ = stroke_plane_offset
        self.air_density_ = 1.18009482370369  # 1.225
        self.wing_length_ = wing_length
        self.mean_chord_ = mean_chord
        self.r33_ = r33
        self.r22_ = r22
        self.r11_ = r11
        self.r00_ = r00
        self.z_cp2_ = z_cp2
        self.z_cp1_ = z_cp1
        self.z_cp0_ = z_cp0
        self.z_rd_ = z_rd
        self.u_ = 0
        self.v_ = 0
        self.w_ = 0
        self.p_ = 0
        self.q_ = 0
        self.r_ = 0
        self.Phi_ = 0
        self.Phi_dot_ = 0
        self.psi_ = 0
        self.psi_dot_ = 0
        self.phi_ = 0
        self.phi_dot_ = 0
        self.theta_ = 0
        self.theta_dot_ = 0
        self.S_Phi_ = 0
        self.C_Phi_ = 0
        self.S_psi_ = 0
        self.C_psi_ = 0
        self.S_phi_ = 0
        self.C_phi_ = 0
        self.U_o1_ = 0
        self.U_o0_ = 0
        self.U_i1_ = 0
        self.U_i0_ = 0
        self.delta_alpha_ = 0
        self.alpha_0_ = 0
        self.alpha_ = 0
        self.C_N_ = 0
        self.a_U2_ = 0
        self.a_U1_ = 0
        self.a_U0_ = 0
        self.K_N2_ = 0
        self.K_N1_ = 0
        self.K_N0_ = 0
        self.K_N_1_ = 0
        self.d_cp_ = 0
        self.r_cp_ = self.R_w_*self.r33_/self.r22_
        self.K_aero2_ = 0
        self.K_aero1_ = 0
        self.K_aero0_ = 0
        self.K_aero_1_ = 0
        self.K_aero_2_ = 0
        self.span_wise_center_of_pressure_ = 0   # 初始化为0
        self.cord_wise_center_of_pressure_ = 0   # 初始化为0
        self.normal_force_ = 0                   # 初始化为0
        self.aero_moment_ = 0                    # 初始化为0
        self.rotational_damping_moment_ = 0      # 初始化为0

    def doNothing(self):
        pass

    def UpdateAeroForce(self):

        # 更新状态量
        self.UpdateVelocityCoeff()
        self.UpdateAoA()
        self.UpdateCN()
        self.UpdateCenterOfPressure()
        self.UpdateVelocitySquaredCoeff()

        # 计算气动力和气动力矩
        self.normal_force_ = 0.5 * self.air_density_ * self.mean_chord_ * self.C_N_ * (
                    self.a_U2_ * pow(self.R_w_, 3) * self.r22_ + self.a_U1_ * pow(self.R_w_,
                                                                                  2) * self.r11_ + self.a_U0_ * self.R_w_ * self.r00_)
        self.aero_moment_ = -0.5 * self.air_density_ * pow(self.mean_chord_, 2) * self.C_N_ * self.d_cp_ * (
                    self.a_U2_ * pow(self.R_w_, 3) * self.z_cp2_ + self.a_U1_ * pow(self.R_w_,
                                                                                    2) * self.z_cp1_ + self.a_U0_ * self.R_w_ * self.z_cp0_)

        # 计算压力中心位置和气动力矩臂
        if self.normal_force_ != 0:
            self.cord_wise_center_of_pressure_ = -self.aero_moment_ / self.normal_force_
        else:
            self.cord_wise_center_of_pressure_ = 0
        self.span_wise_center_of_pressure_ = self.r_cp_

        # 计算旋转阻尼力矩
        self.rotational_damping_moment_ = -0.125 * self.air_density_ * abs(
            self.theta_dot_) * self.theta_dot_ * 5.0 * self.R_w_ * pow(self.mean_chord_, 4) * self.z_rd_

    def UpdateStates(self, body_velocity_roll, body_velocity_pitch, body_velocity_yaw, body_velocity_x, body_velocity_y,
                     body_velocity_z, stroke_plane_angle, stroke_plane_velocity, stroke_angle, stroke_velocity,
                     deviation_angle, deviation_velocity, rotate_angle, rotate_velocity):

        # 更新状态量
        self.u_ = body_velocity_x
        self.v_ = body_velocity_y
        self.w_ = body_velocity_z
        self.p_ = body_velocity_roll
        self.q_ = body_velocity_pitch
        self.r_ = body_velocity_yaw
        self.Phi_ = stroke_plane_angle
        self.Phi_dot_ = stroke_plane_velocity
        self.psi_ = stroke_angle
        self.psi_dot_ = stroke_velocity
        self.phi_ = deviation_angle
        self.phi_dot_ = deviation_velocity
        self.theta_ = rotate_angle
        self.theta_dot_ = rotate_velocity

        # 预先计算三角函数值
        self.S_Phi_ = math.sin(self.Phi_)
        self.C_Phi_ = math.cos(self.Phi_)
        self.S_psi_ = math.sin(self.psi_)
        self.C_psi_ = math.cos(self.psi_)
        self.S_phi_ = math.sin(self.phi_)
        self.C_phi_ = math.cos(self.phi_)

    def UpdateVelocityCoeff(self):
        self.U_o1_ = self.sign_ * (self.p_ * self.C_psi_ * self.C_Phi_ - self.Phi_dot_ * self.S_psi_) + (self.q_ * self.S_psi_ + self.r_ * self.C_psi_ * self.S_Phi_ + self.phi_dot_)
        self.U_o0_ = self.sign_ * (-(self.u_ + self.q_ * self.d_s_) * self.C_phi_ * self.S_Phi_ - self.r_ * self.d_0_ * self.S_phi_ * self.S_psi_ * self.C_Phi_ - (self.v_ - self.p_ * self.d_s_) * self.S_phi_ * self.C_psi_ + self.w_ * self.S_phi_ * self.S_psi_ * self.S_Phi_ + self.p_ * self.d_0_ * self.C_phi_ * self.C_Phi_) + ((self.u_ + self.q_ * self.d_s_) * self.S_phi_ * self.S_psi_ * self.C_Phi_ + self.r_ * self.d_0_ * self.C_phi_ * self.S_Phi_ + self.w_ * self.C_phi_ * self.C_Phi_ + self.p_ * self.d_0_ * self.S_phi_ * self.S_psi_ * self.S_Phi_)
        self.U_i1_ = self.sign_ * (self.p_ * self.S_phi_ * self.S_psi_ * self.C_Phi_ + self.r_ * self.C_phi_ * self.C_Phi_ + self.Phi_dot_ * self.S_phi_ * self.C_psi_) + (-self.p_ * self.C_phi_ * self.S_Phi_ - self.q_ * self.S_phi_ * self.C_psi_ + self.r_ * self.S_phi_ * self.S_psi_ * self.S_Phi_ + self.psi_dot_ * self.C_phi_)
        self.U_i0_ = self.sign_ * (self.r_ * self.d_0_ * self.C_psi_ * self.C_Phi_ - (self.v_ - self.p_ * self.d_s_) * self.S_psi_ - self.w_ * self.C_psi_ * self.S_Phi_) + (-(self.u_ + self.q_ * self.d_s_) * self.C_psi_ * self.C_Phi_ - self.p_ * self.d_0_ * self.C_psi_ * self.S_Phi_)

    def UpdateAoA(self):
        # 计算迎角
        U_i_ = self.U_i1_ * self.r_cp_ + self.U_i0_
        if U_i_ != 0:
            self.delta_alpha_ = math.atan((self.U_o1_ * self.r_cp_ + self.U_o0_) / (self.U_i1_ * self.r_cp_ + self.U_i0_))
        else:
            self.delta_alpha_ = 0
        self.alpha_0_ = self.theta_ + self.sgn(U_i_) * math.pi / 2
        self.alpha_ = self.alpha_0_ - self.delta_alpha_


    def UpdateVelocitySquaredCoeff(self):
        # 更新速度平方系数
        self.a_U2_ = self.U_i1_ * self.U_i1_ + self.U_o1_ * self.U_o1_
        self.a_U1_ = 2 * self.U_i1_ * self.U_i0_ + 2 * self.U_o1_ * self.U_o0_
        self.a_U0_ = self.U_i0_ * self.U_i0_ + self.U_o0_ * self.U_o0_


    def UpdateCN(self):
        self.C_N_ = self.WING_C_N(self.alpha_)

    def UpdateCenterOfPressure(self):
        self.d_cp_ = self.WING_d_cp(self.alpha_)

    def GetSpanCoP(self):
        # 获取横向力矩系数
        return self.span_wise_center_of_pressure_

    def GetChordCoP(self):
        # 获取纵向力矩系数
        return self.cord_wise_center_of_pressure_

    def GetNormalForce(self):
        # 获取法向力系数
        return self.normal_force_

    def GetMoment(self):
        # 获取力矩
        total_moment = self.aero_moment_ + self.rotational_damping_moment_
        return self.rotational_damping_moment_

    def GetM_aero(self):
        # 获取气动力矩
        return self.aero_moment_

    def GetM_rd(self):
        # 获取旋转阻尼力矩
        return self.rotational_damping_moment_

    def WING_C_N(self, alpha):
        # 计算气动力系数
        return 1.8 * math.sin(2 * alpha) * math.cos(alpha) + 1.95 * math.sin(alpha) - 1.5 * math.cos(2 * alpha) * math.sin(alpha)

    def WING_d_cp(self, alpha):
        # 计算压力中心位置
        return 0.46 - 0.332 * math.cos(alpha) - 0.037 * math.cos(3 * alpha) - 0.013 * math.cos(5 * alpha)


    def sgn(self, val):
        return (val > 0) - (val < 0)

    def GetStroke(self):
        return self.psi_

    def GetAoA(self):
        return self.alpha_