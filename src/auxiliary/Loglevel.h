namespace Aux {

class LogLevel {
    public: 
		static void setLevel(int level);
		static int getLevel();
	private:
		static int current;
};

}
