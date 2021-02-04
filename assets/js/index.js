const startBtn = document.querySelector("#start-btn");

const recognition = new webkitSpeechRecognition();
recognition.continuous = true;
recognition.lang = "en-US";
recognition.interimResults = 1;

const synth = window.speechSynthesis;

startBtn.addEventListener("click", () => {
    recognition.start();
});

let utter = new SpeechSynthesisUtterance("Hi, how are you?");
utter.onend = () => {
    recognition.start();
};

recognition.onresult = (e) => {
    const transcript = e.results[e.results.length -1][0].transcript.trim();
    if (transcript === "hello") {
        recognition.stop();
        utter.text = "Hi, how are you?";
        synth.speak(utter);
    }
    else if (transcript === "goodbye") {
        recognition.stop();
        utter.text = "Au Revoir";
        synth.speak(utter);
    }
};