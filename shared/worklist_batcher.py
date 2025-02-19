import copy


class TipState:
    def __init__(self):
        self.tip_state = [0, 0, 0, 0, 0, 0, 0, 0]
        self.tip_contents = [None, None, None, None, None, None, None, None]
        self.max_tip_volume = 400
        self.current_tip_index = 0

    def __str__(self):
        return f"TipState: {self.tip_state}"

    def increment_current_tip_index(self):
        if self.current_tip_index == 7:
            return False
        self.current_tip_index = self.current_tip_index + 1

    def check_if_remaining_tips_are_sufficient(self, volume, content):
        remaining_volume_current_tip = 0
        num_remaining_tips = 0
        if self.tip_contents[self.current_tip_index] == None:
            for i in range(self.current_tip_index, 8):
                num_remaining_tips = num_remaining_tips + 1
        elif self.tip_contents[self.current_tip_index] == content:
            remaining_volume_current_tip = (
                self.max_tip_volume - self.tip_state[self.current_tip_index]
            )
            for i in range(self.current_tip_index + 1, 8):
                num_remaining_tips = num_remaining_tips + 1
        else:
            for i in range(self.current_tip_index + 1, 8):
                num_remaining_tips = num_remaining_tips + 1
        remaining_empty_volume = (
            num_remaining_tips * self.max_tip_volume + remaining_volume_current_tip
        )
        if volume <= remaining_empty_volume:
            return True

        return False

    def add_volume_to_tip(self, volume, content):
        curr_vol = volume
        commands = []
        if self.tip_contents[self.current_tip_index] == None:
            # print("empty", self.current_tip_index)
            for i in range(self.current_tip_index, 8):
                self.tip_state[i] = min(curr_vol, self.max_tip_volume)
                self.tip_contents[i] = content
                commands.append(
                    (self.current_tip_index, min(curr_vol, self.max_tip_volume))
                )
                curr_vol = curr_vol - min(curr_vol, self.max_tip_volume)
                # print("cv", curr_vol)
                if curr_vol - 0.1 <= 0:
                    break
                self.increment_current_tip_index()

        elif self.tip_contents[self.current_tip_index] == content:
            # print("yes", self.current_tip_index)
            init_tip_idx = self.current_tip_index
            for i in range(self.current_tip_index, 8):
                if i == init_tip_idx:
                    v = self.max_tip_volume - self.tip_state[i]
                    self.tip_state[i] = self.tip_state[i] + min(curr_vol, v)
                    commands.append((self.current_tip_index, min(curr_vol, v)))
                    curr_vol = curr_vol - min(curr_vol, v)
                else:
                    # print("hm", content)
                    self.tip_contents[i] = content
                    self.tip_state[i] = self.tip_state[i] + min(
                        curr_vol, self.max_tip_volume
                    )
                    commands.append(
                        (self.current_tip_index, min(curr_vol, self.max_tip_volume))
                    )
                    curr_vol = curr_vol - min(curr_vol, self.max_tip_volume)
                # print("cv", curr_vol)
                if curr_vol - 0.1 <= 0:
                    break
                self.increment_current_tip_index()
        else:
            self.increment_current_tip_index()
            for i in range(self.current_tip_index, 8):
                self.tip_contents[i] = content
                self.tip_state[i] = self.tip_state[i] + min(
                    curr_vol, self.max_tip_volume
                )
                commands.append(
                    (self.current_tip_index, min(curr_vol, self.max_tip_volume))
                )
                curr_vol = curr_vol - min(curr_vol, self.max_tip_volume)
                # print("cv", curr_vol)
                if curr_vol - 0.1 <= 0:
                    break
                self.increment_current_tip_index()

        return commands

    def remove_and_return_remaining_volume_to_tip(self, commands):
        out = []
        tip_idx = 0
        for command in commands:
            vol = command["arg3"]
            while vol > 0:
                if vol - 0.1 <= 0:
                    break
                dis_com = copy.deepcopy(command)
                # print(tip_idx, vol, self.tip_state[tip_idx])
                if vol > self.tip_state[tip_idx]:
                    v_disp = self.tip_state[tip_idx]
                else:
                    v_disp = vol
                dis_com["arg3"] = round(v_disp, 2)
                dis_com["arg4"] = tip_idx
                self.tip_state[tip_idx] = self.tip_state[tip_idx] - v_disp
                vol = vol - v_disp
                # print(dis_com)
                out.append(dis_com)
                if self.tip_state[tip_idx] - 0.1 <= 0:
                    self.tip_contents[tip_idx] = None
                    tip_idx = tip_idx + 1

        return out

    def reset_tip_state(self):
        self.tip_state = [0, 0, 0, 0, 0, 0, 0, 0]
        self.tip_contents = [None, None, None, None, None, None, None, None]
        self.current_tip_index = 0


washer = {
    "action": "wash",
    "arg1": "",
    "arg2": "",
    "arg3": "",
    "arg5": "",
    "arg4": "",
}


def batch_tips(commands):
    tips = TipState()

    all_batches = []
    aspirates = []
    dispenses = []
    batcher = []
    pre_commands = []
    post_commands = []
    for k in commands:
        # print(k)
        if k["action"] == "aspirate":
            loc = k["arg1"] + "_" + k["arg2"]
            if not tips.check_if_remaining_tips_are_sufficient(k["arg3"], loc):
                source_locs = []
                for a in aspirates:
                    source_locs.append(a["arg1"] + "_" + a["arg2"])
                    batcher.append(a)

                # print(dispenses)
                disp_commands = tips.remove_and_return_remaining_volume_to_tip(
                    dispenses
                )
                batcher.extend(disp_commands)
                if len(set(source_locs)) > 1:
                    batcher.append(washer)
                all_batches.append(batcher)
                batcher = []

                aspirates = []
                dispenses = []
                tips.reset_tip_state()

                coms = tips.add_volume_to_tip(k["arg3"], loc)
                for c in coms:
                    if c[1] == 0:
                        continue
                    aspirates.append(
                        {
                            "action": k["action"],
                            "arg1": k["arg1"],
                            "arg2": k["arg2"],
                            "arg3": round(c[1], 2),
                            "arg4": c[0],
                            "arg5": k["arg5"],
                        }
                    )

            else:
                coms = tips.add_volume_to_tip(k["arg3"], loc)
                for c in coms:
                    if c[1] == 0:
                        continue
                    aspirates.append(
                        {
                            "action": k["action"],
                            "arg1": k["arg1"],
                            "arg2": k["arg2"],
                            "arg3": round(c[1], 2),
                            "arg4": c[0],
                            "arg5": k["arg5"],
                        }
                    )

        elif k["action"] == "wash":
            continue
        elif k["action"] == "dispense":
            dispenses.append(k)
        elif k["action"] in ["prepare", "quench", "move"]:
            pre_commands.append(k)
        elif k["action"] in ["seal", "end_quench"]:
            print(k)
            post_commands.append(k)

    if len(aspirates) > 0 or len(dispenses) > 0:
        source_locs = []
        for a in aspirates:
            source_locs.append(a["arg1"] + "_" + a["arg2"])
            batcher.append(a)
        disp_commands = tips.remove_and_return_remaining_volume_to_tip(dispenses)
        batcher.extend(disp_commands)

        if len(set(source_locs)) > 1:
            batcher.append(washer)
        all_batches.append(batcher)

    commands_out = []

    for i in pre_commands:
        commands_out.append(i)

    for i in all_batches:
        for ii in i:
            commands_out.append(ii)

    for i in post_commands:
        commands_out.append(i)

    return commands_out
