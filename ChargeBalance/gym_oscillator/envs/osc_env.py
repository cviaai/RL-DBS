import gym
import numpy as np
from gym import error, spaces, utils
from gym.utils import seeding
import oscillator_cpp
from gym import Env
from gym.spaces import Discrete, MultiDiscrete, MultiBinary, Box


class oscillatorEnv(gym.Env):
  metadata = {'render.modes': ['human']}
  def __init__(self, len_state=250, ep_length=100000, nosc=1000, epsilon=0.03, frrms=0.1,random=False,
               sigmoid=False,charge_balance=False,rectangular=False,init_time=5000,charge_balance_M=0,charge_balance_N=0,integration_step=0.02):

    """ 
    Init function:
    sigmoid: Function that we observe instead of original one: Bool
    len_state: shape of state that agent observes [250,1]: integer

    BVDP params
    nosc: number of oscillators: integer
    epsilon: coupling parameter: float
    frrms: width of the distribution of natural frequencies: float
    """

    super(oscillatorEnv, self).__init__()
  
    self.nosc = nosc
    self.epsilon = epsilon
    self.frrms = frrms
    self.init_time = init_time
    self.ep_length = ep_length
    self.integration_step = integration_step
    #Dimensionality of our observation space
    self.dim = 1
    self.action_space = Box(low=-200, high=200, shape=(1,), dtype=np.float32)
    self.observation_space = Box(low=-1., high=1., shape=(len_state,), dtype=np.float32)

    #Episode Done?
    self.done = False
    self.current_step = 0 
    
    #Our current state, with length(1,len_state)
    self.y_state = []
    self.x_state = []
    #Our actions 
    self.charge_balance = charge_balance
    self.rectangular = rectangular
    #Distance between impulses
    self.charge_balance_M = charge_balance_M
    #Impulse param
    self.charge_balance_N = charge_balance_N
    self.len_state = len_state
    #Reset environment
    self.skip_param = 10
    self.reset()
    
  def step_epsilon(self,action,epsilon):
        #Vectorized form for stable baselines
    if self.epsilon!=epsilon:  
        oscillator_cpp.change_epsilon(epsilon)
        self.epsilon = epsilon
    val = float(action[0])
    oscillator_cpp.Rectangular_signal(val,self.charge_balance_N)
    self.current_step += 1
    
    # if (self.current_step % self.skip_param) == 0:
    self.x_val = oscillator_cpp.Calc_mfx()
    self.y_val = oscillator_cpp.Calc_mfy()

    self.y_state.append(self.y_val)
    self.x_state.append(self.x_val)

    #Check length of our state
    if len(self.y_state) > self.len_state:
        self.y_state = self.y_state[1:]
        self.x_state = self.x_state[1:]

    

    self.done = self.current_step >= self.ep_length

    #Make vectorized form
    arrayed_version = np.array(self.y_state)

    return arrayed_version, self.Reward(self.x_val,self.x_state,val), self.done, {} 
    
  def step(self,action):
        #Vectorized form for stable baselines
    
    val = float(action[0])
    oscillator_cpp.Rectangular_signal(val,self.charge_balance_N)
    self.current_step += 1
    
    # if (self.current_step % self.skip_param) == 0:
    self.x_val = oscillator_cpp.Calc_mfx()
    self.y_val = oscillator_cpp.Calc_mfy()

    self.y_state.append(self.y_val)
    self.x_state.append(self.x_val)

    #Check length of our state
    if len(self.y_state) > self.len_state:
        self.y_state = self.y_state[1:]
        self.x_state = self.x_state[1:]

    

    self.done = self.current_step >= self.ep_length

    #Make vectorized form
    arrayed_version = np.array(self.y_state)

    return arrayed_version, self.Reward(self.x_val,self.x_state,val), self.done, {} 





  def pertrubation(self,action):
        #Vectorized form for stable baselines
    
    val = float(action[0])
    oscillator_cpp.Pertrubation(val)
    self.current_step += 1
    
    # if (self.current_step % self.skip_param) == 0:
    self.x_val = oscillator_cpp.Calc_mfx()
    self.y_val = oscillator_cpp.Calc_mfy()

    self.y_state.append(self.y_val)
    self.x_state.append(self.x_val)
  
    #Check length of our state
    if len(self.y_state) > self.len_state:
        self.y_state = self.y_state[1:]
        self.x_state = self.x_state[1:]

    

    self.done = self.current_step >= self.ep_length

    #Make vectorized form
    arrayed_version = np.array(self.y_state)

    return arrayed_version, self.Reward(self.x_val,self.x_state,val), self.done, {} 

  def reset(self):
    """
    Reset environment, and get a window 250 of self.len_state size

    Returns:arrayed_version:np.array(1,len_state)

    """
    self.current_step = 0 
    self.y_state = []
    self.x_state = []
    self.actions = []
    oscillator_cpp.init(self.nosc,self.epsilon,self.frrms,self.init_time,np.random.randint(1,10000),self.integration_step)

    for i in range(self.len_state):
        oscillator_cpp.Make_step()
        
        if (self.current_step % self.skip_param) == 0:
            self.x_val = oscillator_cpp.Calc_mfx()
            self.y_val = oscillator_cpp.Calc_mfy()

            self.y_state.append(self.y_val)
            self.x_state.append(self.x_val)
            
        #self.x_val = oscillator_cpp.Calc_mfx()
        #self.y_val = oscillator_cpp.Calc_mfy()
        
        #self.y_state.append(self.y_val)
        #self.x_state.append(self.x_val)
        #self.actions.append(0)

        #Check length of our state
        if len(self.y_state) > self.len_state:
            self.y_state = self.y_state[1:]
            self.x_state = self.x_state[1:]

    arrayed_version = np.array(self.y_state)    
    
    return arrayed_version
    
  def render(self, mode='human', close=False):
    """
    Pass...
    """
    pass

  def Reward(self, x,x_state,action_value,baseline = False):
    """
    Super duper reward function, i am joking, just sum of absolute values which we supress + penalty for actions
    returns: float
    """
    
    return -(x-np.mean(x_state))**2 - (1e-2)*np.abs(action_value)


