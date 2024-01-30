#include "controller.h"

CController::CController()
{
	_k = 7;
	Initialize();
	j = 0;
}

CController::~CController()
{
}



void CController::read(double t, double* q, double* qdot)
{	
	_t = t;
	if (_bool_init == true)
	{
		_init_t = _t;
		_bool_init = false;
	}

	_dt = t - _pre_t;

	_pre_t = t;

	for (int i = 0; i < _k; i++)
	{
		_q(i) = q[i];
		_qdot(i) = qdot[i];
		_pre_q(i) = _q(i);
		_pre_qdot(i) = _qdot(i);		
		if(_t < 2.0)
        {
			_qdot(i) = qdot[i];
		}
	}
}

void CController::write(double* torque)
{
	for (int i = 0; i < _k; i++)
	{
		torque[i] = _torque(i);
	}
}

void CController::control_mujoco()
{
    ModelUpdate();
    motionPlan();
	
	if(_control_mode == 1) // Joint Space
	{
		if (_t - _init_t < 0.1 && _bool_joint_motion == false)
		{
			_start_time = _init_t;
			_end_time = _start_time + _motion_time;
			JointTrajectory.reset_initial(_start_time, _q, _qdot);
			JointTrajectory.update_goal(_q_goal, _qdot_goal, _end_time);
			_bool_joint_motion = true;
		}
		
		JointTrajectory.update_time(_t);
		_q_des = JointTrajectory.position_cubicSpline();
		_qdot_des = JointTrajectory.velocity_cubicSpline();

		JointControl();

		if (JointTrajectory.check_trajectory_complete() == 1)
		{
			_bool_plan(_cnt_plan) = 1;
			j++;
			cout << "joint space complete " << j<<endl;
			cout<<"jointspace _x_hand 	: "<<_x_hand.transpose()<<endl;
			_bool_init = true;
		}
	}
	else if(_control_mode == 2) // Operational Space
	{		
		if (_t - _init_t < 0.1 && _bool_ee_motion == false)
		{
			_start_time = _init_t;
			_end_time = _start_time + _motion_time;
			HandTrajectory.reset_initial(_start_time, _x_hand, _xdot_hand);
			HandTrajectory.update_goal(_x_goal_hand, _xdot_goal_hand, _end_time);
			_bool_ee_motion = true;
			// cout<<"_t : "<<_t << _init_t <<endl;
			// cout<<"_x_hand 	: "<<_x_hand.transpose()<<endl;
		}

		HandTrajectory.update_time(_t);
		_x_des_hand.head(3) = HandTrajectory.position_cubicSpline();
		_R_des_hand = HandTrajectory.rotationCubic();
		_x_des_hand.segment<3>(3) = CustomMath::GetBodyRotationAngle(_R_des_hand); // rpy
		_xdot_des_hand.head(3) = HandTrajectory.velocity_cubicSpline();
		_xdot_des_hand.segment<3>(3) = HandTrajectory.rotationCubicDot();		
		
		if (HandTrajectory.check_trajectory_complete() == 1)
		{
			
			_bool_plan(_cnt_plan) = 1;
			_bool_init = true;
			j++;
			cout << "OSC Complete " << j << endl;
			cout<<"OSC _x_hand 	: "<<_x_hand.transpose()<<endl;
		}
		OperationalControl();
		
	}
		
}

void CController::ModelUpdate()
{
    Model.update_kinematics(_q, _qdot);
	Model.update_dynamics();
    Model.calculate_EE_Jacobians();
	Model.calculate_EE_positions_orientations();
	Model.calculate_EE_velocity();

	_J_hands = Model._J_hand;
	_Jdot_qdot = _J_hands * _qdot;

	_x_hand.head(3) = Model._x_hand;
	_x_hand.tail(3) = CustomMath::GetBodyRotationAngle(Model._R_hand);
	_xdot_hand = Model._xdot_hand;

}	

void CController::motionPlan()
{	
	if (_bool_plan(_cnt_plan) == 1)
	{

		if(_cnt_plan == 0)
		{
			cout << "cnt_plan = 0" << endl; 
			_q_order(0) = -1.47;
			_q_order(1) = -0.5;
			_q_order(2) = -0.4;
			_q_order(3) = -1.3;
			_q_order(4) = 0.3;
			_q_order(5) = 1.0;
			_q_order(6) = -0.1;		                    
			reset_target(7.0, _q_order, _qdot);
			_cnt_plan++;
		}
		else if(_cnt_plan == 1)
		{
			cout << "cnt_plan = 1" << endl; 
			goal_pos(0) = 0.2;
			goal_pos(1) = 0.0;
			goal_pos(2) = 0.78;
			goal_ori(0) = -2.4;
			goal_ori(1) = 0;
			goal_ori(2) = -1.0;

			reset_target(10.0, goal_pos, goal_ori);
			_cnt_plan++;
		}
		
	}
}

// Task Space
void CController::reset_target(double motion_time, Vector3d target_pos, Vector3d target_ori)
{
	_control_mode = 2;
	_motion_time = motion_time;
	_bool_joint_motion = false;
	_bool_ee_motion = false;

	_x_goal_hand.head(3) = target_pos;
	_x_goal_hand.tail(3) = target_ori;
	_xdot_goal_hand.setZero();
}

// Joint Space
void CController::reset_target(double motion_time, VectorXd target_joint_position, VectorXd target_joint_velocity)
{
	_control_mode = 1;
	_motion_time = motion_time;
	_bool_joint_motion = false;
	_bool_ee_motion = false;

	_q_goal = target_joint_position.head(7);
	// _qdot_goal = target_joint_velocity.head(7);
	_qdot_goal.setZero();
}

void CController::JointControl()
{	
	_torque.setZero();
	_A_diagonal = Model._A;
	for(int i = 0; i < 7; i++){
		_A_diagonal(i,i) += 1.0;
	}
	_torque = _A_diagonal*(400*(_q_des - _q) + 40*(_qdot_des - _qdot)) + Model._bg;
}

void CController::CLIK()
{
	_torque.setZero();	
	_dt = 0.003;

	_x_err_hand.segment(0,3) = _x_des_hand.head(3) - _x_hand.head(3);
	_x_err_hand.segment(3,3) = -CustomMath::getPhi(Model._R_hand, _R_des_hand);

	_J_bar_hands = CustomMath::pseudoInverseQR(_J_hands);

	_qdot_des = _J_bar_hands*(_xdot_des_hand + _x_kp*(_x_err_hand));
	_q_des = _q_des + _dt*_qdot_des;
	_A_diagonal = Model._A;
	for(int i = 0; i < 7; i++){
		_A_diagonal(i,i) += 1.0;
	}

	_torque = _A_diagonal * (400 * (_q_des - _q) + 40 * (_qdot_des - _qdot)) + Model._bg;

}


void CController::OperationalControl()
{
	_torque.setZero();
	_J_T_hands.setZero(_k, 6);

	_x_err_hand.segment(0,3) = _x_des_hand.head(3) - _x_hand.head(3); // xyz
	_x_err_hand.segment(3,3) = -CustomMath::getPhi(Model._R_hand, _R_des_hand); // orientation
	_xdot_err_hand.segment(0,3) = _xdot_des_hand.head(3) - _xdot_hand.head(3); // xyz
	_xdot_err_hand.segment(3,3) = _xdot_des_hand.segment(3,3) - _xdot_hand.segment(3,3); // orientation

	_J_T_hands = _J_hands.transpose(); // 7x6
	_J_bar_hands = CustomMath::pseudoInverseQR(_J_hands); // weighted pseudo inverse
	_J_barT_hands = _J_bar_hands.transpose(); // 6x7
	_Null = _I - _J_T_hands * _J_barT_hands; // 7x7
	

	_lambda = _J_barT_hands * Model._A * _J_bar_hands; // 6x6
	_F_star = _kpj * _x_err_hand + _kdj * _xdot_err_hand; // 6x1
	
	_q_des0 = _q;
	_q_des0(0) = -70 * DEG2RAD;
	_second_task = 30 * (_q_des0 - _q);


	_torque = _J_T_hands * _lambda * _F_star + Model._bg +  _Null *Model._A * _second_task; // 7x1
	// _torque = _J_T_hands * _lambda * _F_star + Model._bg; // 7x1

	// cout << "torque:" << _torque << "\n";
}


void CController::Initialize()
{
	
	j = 0;
    _control_mode = 1; //1: joint space, 2: task space(CLIK)

	_bool_init = true;
	_t = 0.0;
	_init_t = 0.0;
	_pre_t = 0.0;
	_dt = 0.0;
	plan_b = false;
	_kpj = 400.0;
	_kdj = 20.0;

	_kps = 70;

	_x_kp =1;//작게 0.1

    _q.setZero(_k);
	_qdot.setZero(_k);
	_torque.setZero(_k);

	_J_hands.setZero(6,_k);
	_J_bar_hands.setZero(_k,6);

	_x_hand.setZero(6);
	_xdot_hand.setZero(6);

	_bool_plan.setZero(30);

	_q_home.setZero(_k);
	_q_home(0) = 0.374;
	_q_home(1) = -1.02;
	_q_home(2) = 0.245;
	_q_home(3) = -1.51;
	_q_home(4) = 0.0102;
	_q_home(5) = 0.655;
	_q_home(6) = 0.3;

	_start_time = 0.0;
	_end_time = 0.0;
	_motion_time = 0.0;

	_bool_joint_motion = false;
	_bool_ee_motion = false;

	_q_des.setZero(_k);
	_qdot_des.setZero(_k);
	_q_goal.setZero(_k);
	_qdot_goal.setZero(_k);

	_x_des_hand.setZero(6);
	_xdot_des_hand.setZero(6);
	_x_goal_hand.setZero(6);
	_xdot_goal_hand.setZero(6);
	_Jdot_qdot.setZero(6);

	_pos_goal_hand.setZero(); // 3x1 
	_rpy_goal_hand.setZero(); // 3x1
	JointTrajectory.set_size(_k);
	_A_diagonal.setZero(_k,_k);

	_x_err_hand.setZero(6);
	_R_des_hand.setZero();

	_I.setIdentity(7,7);

	_pre_q.setZero(7);
	_pre_qdot.setZero(7);

	_q_order.setZero(7);
	_qdot_order.setZero(7);
	_cnt_plan = 0;
	_bool_plan(_cnt_plan) = 1;


	///////////OperationalSpaceControl//////////////
	goal_pos.setZero(3,1);
	goal_ori.setZero(3,1);
	goal_.setZero(6,1);
	_J_T_hands.setZero(_k, 6);
	_J_barT_hands.setZero(6, _k);
	_Null.setZero(_k, _k);
	_lambda.setZero(6, 6);
	_x_err_hand.setZero(6);
	_xdot_err_hand.setZero(6);

	_F_star.setZero(_k);
	_q_des0.setZero(_k);
	_second_task.setZero(_k);
}
