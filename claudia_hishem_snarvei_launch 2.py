def verify():

    from ast2000tools import constants as con
    from ast2000tools import solar_system
    from ast2000tools import space_mission
    from ast2000tools import shortcuts
    import numpy as np
    from ast2000tools import utils
    utils.check_for_newer_version()

    print( utils.get_seed('Claudieg') )
    seed = 93823# insert student's seed number
    # sjekk at dette er hva dere får ^

    code_escape_trajectory =  45499
    code_launch_results    =  69502

    """------------------------------------------------------------------"""

    mission = space_mission.SpaceMission(seed)
    system = solar_system.SolarSystem(seed)

    shortcut = shortcuts.SpaceMissionShortcuts(mission, [code_escape_trajectory,code_launch_results])

    ################################################################
    #            PLACE SPACECRAFT ON ESCAPE TRAJECTORY             #
    ################################################################
    #                   |      For Part 3      |
    #                   ------------------------

    # OBS!!! Here the angle is in degrees!
    # Copy and paste from here, but choose only one of the two parts (A or
    # B) before copying the rest. Parts A and B depends on what the
    # student(s) actually need this shortcut for. The parts are:

    # SHORTCUT A:
    # use this if the student(s) only need shortcut for generalised launch,
    # i.e., were able to launch their spacecraft in Part 1.

    # SHORTCUT B:
    # use this if the student(s) needed Shortcut B in the two previous
    # shortcuts.
    """
    ------------------------------------------------------------------------
    """




    """
    Documentation
    place_spacecraft_on_escape_trajectory(
        rocket_thrust, rocket_mass_loss_rate, time height_above_surface,
        direction_angle, remaining_fuel_mass):

    ------------------------------------------------------------------------
    place_spacecraft_on_escape_trajectory() places the spacecraft on an
    escape trajectory pointing directly away from the home planet.

    Parameters
    ----------
    rocket_thrust  :  float
        The total thrust of the rocket, in NEWTONS.

    rocket_mass_loss_rate  :  float
        The total mass loss rate of the rocket, in KILOGRAMS PER SECOND.

    time  :  float
        The time at which the spacecraft should be placed on the escape
        trajectory, in YEARS from the initial system time.

    height_above_surface  :  float
        The heigh above the home planet surface to place the spacecraft, in
        METERS (after launch).

    direction_angle  :  float
        The angle of the direction of motion of the spacecraft with respect
        to the x-axis, in DEGREES.

    remaining_fuel_mass  :  float
        The mass of fuel carried by the spacecraft after placing it on the
        escape trajectory, in KILOGRAMS.

    Raises
    ------
    RuntimeError
        When none of the provided codes are valid for unlocking this method.
    ------------------------------------------------------------------------

    """

    # ----------
    # Shortcut A
    # ----------

    thrust = 1500000 # insert the thrust force of your spacecraft here
    mass_loss_rate = 200 # insert the mass loss rate of your spacecraft here

    # ----------
    # Shortcut B (if the student(s) need compute_engine_performance())
    # ----------


    """
    GROUP TEACHER OF AST2000:
    ------------------------------------------------------------------------
    """
    # Use shortcut B of the shortcut compute_engine_performance for thrust
    # and mass_loss_rate. The students that need this part should already
    # have this when they used the shortcut for Part 1.B. 
    """
    ------------------------------------------------------------------------
    """


    # choose these values freely, but they should be relevant to where you
    # want to go, e.g., if you want to travel outwards of your solar system,
    # don't let the direction angle be 0 if you are launching from
    # coordinates close to (-x, 0), as this will send you in the opposite
    # direction), and vice versa if your destination is a planet closer to
    # your sun
    '''Timestep from hohmann transfer function "Del5_1"
    Angle launches about same direction as planet trajectory at timestep
    Height above surface is approx where we are after escape velocity is reached'''
    time =  76965*1e-4# insert the time for the spacecraft to be put on escape trajectory
    height_above_suface = 10386466 # insert height above surface you want the rocket placed
    direction_angle = 320# insert the angle between the x-axis and the rocket's motion
    fuel_left = 30000# insert how much fuel you want for your trip

    shortcut.place_spacecraft_on_escape_trajectory(
        thrust,
        mass_loss_rate,
        time,
        height_above_suface,
        direction_angle,
        fuel_left
        )

    """
    GROUP TEACHER OF AST2000:
    ------------------------------------------------------------------------
    """
    # In order to get the position of the spacecraft after launch, the
    # students will have to use this shortcut from Part 1 of the project,
    # i.e., this should not affect their score for Part 3.
    """
    ------------------------------------------------------------------------
    """

    fuel_consumed, time_after_launch, pos_after_launch, vel_after_launch\
        = shortcut.get_launch_results()
    # fuel_consumed is None because place_spacecraft_on_escape_trajectory()
    # don't actually launch the rocket, but just place the rocket in correct
    # position and with the correct velocity.
    mission.verify_launch_result(pos_after_launch)
    mission.take_picture()
    # print(mission.measure_star_doppler_shifts())
    # mission.verify_manual_orientation(pos_after_launch,vel_after_launch,199)
    # print(pos_after_launch,vel_after_launch)
    print(pos_after_launch)
    print(vel_after_launch)
    # print(time_after_launch)
    # print(fuel_consumed)
if __name__ == '__main__':

    verify()
