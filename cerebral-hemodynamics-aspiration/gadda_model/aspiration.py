"""
Venous aspiration protocols

Closed-loop reinfusion system with pressure-driven flow.
Includes constant-flow, stepwise escalation, and dose-response protocols.
"""

from typing import List, Literal
import numpy as np


class PIDController:
    """
    PID controller for ICP regulation during venous aspiration.
    
    Implements proportional-integral-derivative control to maintain target ICP
    by adjusting aspiration flow rate based on ICP error.
    """
    
    def __init__(self, target_ICP: float, Kp: float = 0.5, Ki: float = 0.01, 
                 Kd: float = 0.05, dt: float = 1.0, min_flow: float = 0, 
                 max_flow: float = 150):
        """
        Initialize PID controller.
        
        Args:
            target_ICP: Target intracranial pressure (mmHg)
            Kp: Proportional gain
            Ki: Integral gain
            Kd: Derivative gain
            dt: Time step (seconds)
            min_flow: Minimum flow rate (mL/min)
            max_flow: Maximum flow rate (mL/min)
        """
        self.target_ICP = target_ICP
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.dt = dt
        self.min_flow = min_flow
        self.max_flow = max_flow
        self.integral = 0.0
        self.prev_error = 0.0
        self.last_output = 0.0  # For anti-windup
    
    def __call__(self, t: float, ICP_current: float) -> float:
        """
        Calculate control output (flow rate in mL/s).
        
        Args:
            t: Current time (seconds)
            ICP_current: Current ICP measurement (mmHg)
            
        Returns:
            Flow rate in mL/s
        """
        # Error = Current - Target (positive error means ICP too high, need more flow)
        error = ICP_current - self.target_ICP
        

        self.integral += error * self.dt
        # Clamp output to prevent integral windup
        self.integral = np.clip(self.integral, -100.0, 100.0)
        
        # Calculate derivative
        derivative = (error - self.prev_error) / self.dt if self.dt > 0 else 0.0
        
        # PID output
        flow_mL_per_min = self.Kp * error + self.Ki * self.integral + self.Kd * derivative
        
        # Clip to safety limits
        flow_mL_per_min = np.clip(flow_mL_per_min, self.min_flow, self.max_flow)
        

        self.prev_error = error
        
        return flow_mL_per_min / 60.0  # Convert to mL/s


class AspirationProtocol:
    """Aspiration protocol manager for venous aspiration interventions."""
    
    def __init__(self, protocol_dict):
        self.protocol = protocol_dict
    
    def __getitem__(self, key):
        return self.protocol[key]
    
    def __setitem__(self, key, value):
        self.protocol[key] = value
    
    def get(self, key, default=None):
        return self.protocol.get(key, default)
    
    def get_flow_rate(self, t):
        """Get flow rate at time t (mL/s)."""
        for segment in self.protocol['flow_schedule']:
            if segment['start'] <= t < segment['end']:
                frac = (t - segment['start']) / (segment['end'] - segment['start'])
                flow = segment['flow_start'] + frac * (segment['flow_end'] - segment['flow_start'])
                return flow
        return 0.0
    
    @staticmethod
    def create_constant_flow(
        site: Literal['Pvs', 'Pv', 'J3', 'J2', 'J1'],
        target_flow_mL_min: float = 120.0,
        ramp_duration_s: float = 30.0,
        total_duration_s: float = 600.0
    ) -> 'AspirationProtocol':
        """
        Create constant flow aspiration protocol.
        
        Args:
            site: Aspiration location ('Pvs', 'Pv', 'J3', 'J2', 'J1')
            target_flow_mL_min: Target flow rate (mL/min), default 120
            ramp_duration_s: Ramp-up time (seconds), default 30
            total_duration_s: Total protocol duration (seconds), default 600 (10 min)
        
        Returns:
            AspirationProtocol instance
        """
        target_flow_mL_s = target_flow_mL_min / 60.0
        
        protocol_dict = {
            'site': site,
            'mode': 'piecewise',
            'flow_schedule': [
                {
                    'start': 0.0,
                    'end': ramp_duration_s,
                    'flow_start': 0.0,
                    'flow_end': target_flow_mL_s,
                    'ramp': 'linear'
                },
                {
                    'start': ramp_duration_s,
                    'end': total_duration_s,
                    'flow_start': target_flow_mL_s,
                    'flow_end': target_flow_mL_s,
                    'ramp': 'step'
                }
            ]
        }
        
        return AspirationProtocol(protocol_dict)
    
    @staticmethod
    def create_stepwise_escalation(
        site: Literal['Pvs', 'Pv', 'J3', 'J2', 'J1'],
        flow_steps_mL_min: List[float] = None,
        step_duration_s: float = 180.0,
        ramp_between_steps_s: float = 30.0
    ) -> 'AspirationProtocol':
        """
        Create stepwise escalation protocol.
        
        Args:
            site: Aspiration location
            flow_steps_mL_min: List of flow rates (mL/min) for each step
                Default: [60, 120, 180, 240]
            step_duration_s: Duration of each step (seconds), default 180 (3 min)
            ramp_between_steps_s: Ramp time between steps (seconds), default 30
        
        Returns:
            AspirationProtocol instance
        """
        if flow_steps_mL_min is None:
            flow_steps_mL_min = [60, 120, 180, 240]
        
        flow_schedule = []
        t_current = 0.0
        
        # Initial ramp to first step
        flow_0_mL_s = flow_steps_mL_min[0] / 60.0
        flow_schedule.append({
            'start': 0.0,
            'end': ramp_between_steps_s,
            'flow_start': 0.0,
            'flow_end': flow_0_mL_s,
            'ramp': 'linear'
        })
        t_current = ramp_between_steps_s
        
        # Each step
        for i, flow_mL_min in enumerate(flow_steps_mL_min):
            flow_mL_s = flow_mL_min / 60.0
            

            flow_schedule.append({
                'start': t_current,
                'end': t_current + step_duration_s,
                'flow_start': flow_mL_s,
                'flow_end': flow_mL_s,
                'ramp': 'step'
            })
            t_current += step_duration_s
            
            # Ramp to next step (if not last)
            if i < len(flow_steps_mL_min) - 1:
                flow_next_mL_s = flow_steps_mL_min[i + 1] / 60.0
                flow_schedule.append({
                    'start': t_current,
                    'end': t_current + ramp_between_steps_s,
                    'flow_start': flow_mL_s,
                    'flow_end': flow_next_mL_s,
                    'ramp': 'linear'
                })
                t_current += ramp_between_steps_s
        
        protocol_dict = {
            'site': site,
            'mode': 'piecewise',
            'flow_schedule': flow_schedule
        }
        
        return AspirationProtocol(protocol_dict)


__all__ = ['AspirationProtocol', 'PIDController']
