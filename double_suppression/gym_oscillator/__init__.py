from gym.envs.registration import register

register(
    id='oscillator-v0',
    entry_point='gym_oscillator.envs:oscillatorEnv',
)
