template<class Derived>
class state_machine
{

protected:
	template<int CurrentState, class Event, int NextState, void (Derived::*action)(Event const&)>
	struct transition
	{
		// for later use by our metaprogram
		static int const current_state = CurrentState;
		static int const next_state = NextState;
		typedef Event event;
		typedef Derived fsm_t;
		// do the transition action
		static void execute(Derived& fsm, Event const& e)
		{
			(fsm.*action) (e);
		}
	};
};


// concrete FSM implementation
class player : public state_machine< player>
{
private:
	// the list of FSM states
	enum states { Empty, Open, Stopped, Playing, Paused , initial_state = Empty };
	
	// transition actions
	void start_playback(play const&);
	void open_drawer(open_close const&);
	void close_drawer(open_close const&);
	void store_cd_info(cd_detected const&);
	void stop_playback(stop const&);
	void pause_playback(pause const&);
	void resume_playback(play const&);
	void stop_and_open(open_close const&);
	
	friend class state_machine<player>;

	// transition table
	typedef player p; // makes transition table cleaner
	struct transition_table : mpl::vector11<
		transition < Stopped , play 	, Playing , &p::start_playback >,
		transition < Stopped , open_close , Open , &p::open_drawer >,
		// +---------+-------------+---------+---------------------+
		transition < Open , open_close , Empty , &p::close_drawer >,
		// +---------+-------------+---------+---------------------+
		transition < Empty , open_close , Open , &p::open_drawer >,
		transition < Empty , cd_detected , Stopped , &p::store_cd_info >,
		// +---------+-------------+---------+---------------------+
		transition < Playing , stop 	, Stopped , &p::stop_playback >,
		transition < Playing , pause 	, Paused , &p::pause_playback >,
		transition < Playing , open_close , Open , &p::stop_and_open >,
		// +---------+-------------+---------+---------------------+
		transition < Paused , play 	, Playing , &p::resume_playback >,
		transition < Paused , stop 	, Stopped , &p::stop_playback >,
		transition < Paused , open_close , Open , &p::stop_and_open >
		// +---------+-------------+---------+---------------------+
	> {};
};

