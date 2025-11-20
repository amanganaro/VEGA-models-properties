package insilico.dio1.model.modelGBTpmml;

import java.util.ArrayList;
import java.util.HashMap;

public class XmlTagUtils {

    public static ArrayList<String> divideTags(String xml) {
        ArrayList<String> res = new ArrayList<>();

        int idx = xml.toLowerCase().indexOf("<");
        if (idx == -1)
            return res;

        String curTag = "";
        while (idx < xml.length()) {
            if (xml.charAt(idx) == '<') {
                curTag = "" + xml.charAt(idx);
            } else if (xml.charAt(idx) == '>') {
                curTag += xml.charAt(idx);
                res.add(curTag);
            } else {
                curTag += xml.charAt(idx);
            }
            idx++;
        }

        return res;
    }


    public static String getTagName(String xml) {
        String res = "";
        int idx = xml.toLowerCase().indexOf("<");
        if (idx == -1)
            return res;
        idx++;

        while (idx < xml.length()) {
            if ( (xml.charAt(idx) == ' ') || (xml.charAt(idx) == '>') )
                return res;
            else
                res += xml.charAt(idx);
            idx++;
        }

        return res;
    }


    public static boolean isTagClosed(String xml) {
        int idx = xml.toLowerCase().indexOf(">");
        if (idx == -1)
            return false;
        if (xml.charAt(idx-1) == '/')
            return true;
        return false;
    }


    public static String composeTags(ArrayList<String> tags) {
        String res = "";
        for (String s : tags)
            res += s;
        return res;
    }

    public static HashMap<String, String> getClosedTag(String xml, String tagName) {

        ArrayList<String> tags = divideTags(xml);

        ArrayList<String> res = new ArrayList<>();
        ArrayList<String> remaining = new ArrayList<>();

        int nOpen = 0;
        boolean begin = false;
        boolean ended = false;
        for (String curTag : tags) {

            if (ended)
                remaining.add(curTag);

            if (getTagName(curTag).equalsIgnoreCase(tagName)) {
                begin = true;
                res.add(curTag);
                if (!isTagClosed(curTag))
                    nOpen++;
                else if (nOpen == 0) {
                    ended = true;
                }
                continue;
            }

            if (!begin)
                continue;

            if (getTagName(curTag).equalsIgnoreCase("/" + tagName)) {
                res.add(curTag);
                nOpen--;
                if (nOpen == 0)
                    ended = true;
                continue;
            }

            res.add(curTag);
        }

        HashMap<String, String> r = new HashMap<>();
        r.put("result", composeTags(res));
        r.put("remaining", composeTags(remaining));
        return r;
    }

}
