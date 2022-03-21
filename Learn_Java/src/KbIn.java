public class KbIn {
    public static void main(String[] args) 
        throws java.io.IOException{ // this line handles input error
            char ch, ignore, answer = 'K';

            do{
                System.out.println("I'm thinking of a letter between A and Z.");
                System.out.print("Try to guess it: ");

                ch = (char) System.in.read(); // read a char

                // discard any other char in the input buffer
                // Pressing ENTER causes a carriage return and a line feed (newline) sequence to be generated.
                do {
                    ignore = (char) System.in.read();
                } while (ignore != '\n'); 

                if (ch == answer) System.out.println("** Right! **");
                else { 
                    System.out.println("Sorry, you're wrong :(");
                    if ( ch < answer) System.out.println("too low");
                    else System.out.println("too high");
                    System.out.println("Try again!\n");
                }
            } while (answer != ch);
    }    
}
